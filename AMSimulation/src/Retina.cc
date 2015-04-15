#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Retina.h"
using namespace slhcl1tt;

#include <iostream>
#include <string>

#include "TMath.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"

Retina::Retina(std::vector<xyPoint> hits_, unsigned int pbins_, unsigned int qbins_,
               double pmin_, double pmax_, double qmin_, double qmax_,
               std::vector<double> sigma_, double minWeight_, FitView view_) :
  hits(hits_),
  pbins(pbins_),
  qbins(qbins_),
  pmin(pmin_),
  pmax(pmax_),
  offset(0.),
  qmin(qmin_+offset),
  qmax(qmax_+offset),
  sigma(sigma_),
  minWeight(minWeight_),
  view(view_)
{

  pbinsize = (pmax-pmin)/(double)pbins;
  qbinsize = (qmax-qmin)/(double)qbins;

  makeGrid();
}

Retina::~Retina() {
}


void Retina::makeGrid() {
  Grid.resize(pbins);
  for (unsigned int i = 0; i < pbins; i++) {
    Grid[i].resize(qbins);
  }
}

void Retina::fillGrid() {
  for (unsigned int i = 0; i < pbins; i++) {
    for (unsigned int j = 0; j < qbins; j++) {
      double p_value = pmin + pbinsize * (i + 0.5);
      double q_value = qmin + qbinsize * (j + 0.5);
      Grid[i][j] = getResponseXpXm(p_value, q_value);
      //std::cout << i << " " << j << " p = " << p_value << " q = " << q_value << " w = " << Grid[i][j] << std::endl;
    }
  }
}

double Retina::getResponseXpXm(double x_plus, double x_minus) {

  double Rij = 0.;
  unsigned int hits_tot = hits.size();

  // *** This is always assuming barrel hits, it has to be checked
  double y0 = 23.0;
  double y1 = 108.0;
  if ( view == XY ){
    y0 = 0.5/y0;
    y1 = 0.5/y1;
  }
  // ***

  // Get p and q from x_plus and x_minus
  const double x0 = x_plus - x_minus;
  const double x1 = x_plus + x_minus;
  const double p = (y1 - y0) / (x1 - x0);
  const double q = y0 - p*x0;

  for (unsigned int kr = 0; kr < hits_tot; kr++) {

    double x = hits[kr].x;
    double y = hits[kr].y;

    double sigma_local = sigma[0];

    if ( view == RZ && fabs(y) > 60. )
      sigma_local = sigma[3];

    //double Sijkr = (p*x-y+q)/p;  // from  x - (y-q)/p
    double Sijkr = fabs(y-p*x-q)/sqrt(1.+p*p);

    double term = exp(-0.5*Sijkr*Sijkr/(sigma_local*sigma_local));

    Rij += term;
  }
  if (Rij < 1e-6) return 1e-6;
  else return Rij;

}


void Retina::dumpGrid(int eventNum, int step, int imax) {

//  std::cout << std::endl;
//  std::cout << "Dumping the retina ... " << std::endl;
//  std::cout << " p range = " << pmin << " --> " << pmax << std::endl;
//  std::cout << "   p bins = " << pbins << std::endl;
//  std::cout << "   p step = " << pbinsize << std::endl;
//  std::cout << " q range = " << qmin << " --> " << qmax << std::endl;
//  std::cout << "   q bins = " << qbins << std::endl;
//  std::cout << "   q step = " << qbinsize << std::endl;
//  std::cout << " sigma = " << sigma[0] << std::endl;
//  std::cout << " minWeight = " << minWeight << std::endl;
//  for (unsigned int i = 0; i < pbins; i++){
//    for (unsigned int j = 0; j < qbins; j++)
//      std::cout <<  Grid[i][j] << "\t";
//    std::cout << std::endl;
//  }
//  std::cout << std::endl;


  TH2D pq_h("pq_h", "x_{+}-x_{-} grid;x_{+};x_{-}", pbins, pmin, pmax, qbins, qmin, qmax);
  for (unsigned int i = 0; i < pbins; i++) {
    for (unsigned int j = 0; j < qbins; j++) {
      pq_h.SetBinContent(i+1, j+1, Grid[i][j]);
    }
  }
  gStyle->SetPalette(53);
  gStyle->SetPaintTextFormat("5.2f");
  TCanvas c("c", "c", 650, 600);
  c.SetRightMargin(0.1346749);
  pq_h.SetStats(false);
  pq_h.SetMaximum(7.);
  pq_h.SetMinimum(0.);
  std::string draw_s("COLZTEXT");
  if (qbins > 30 || pbins > 30) draw_s = "COLZ";
  pq_h.Draw(draw_s.c_str());
  std::string fit_view = "XY";
  if (view == RZ) fit_view = "RZ";
  c.SaveAs(Form("PQgrid_%d_%s_%d-%d.png",eventNum,fit_view.c_str(),step,imax));

}


void Retina::findMaxima() {

  for (unsigned int i = 1; i < pbins-1; i++) {
    for (unsigned int j = 1; j < qbins-1; j++) {
      if (    Grid[i][j] > Grid[i-1][j]
           && Grid[i][j] > Grid[i+1][j]
           && Grid[i][j] > Grid[i][j-1]
           && Grid[i][j] > Grid[i][j+1]

           && Grid[i][j] > Grid[i+1][j+1]
           && Grid[i][j] > Grid[i+1][j-1]
           && Grid[i][j] > Grid[i-1][j-1]
           && Grid[i][j] > Grid[i-1][j+1]
         ) {
        if (Grid[i][j] < minWeight) continue; // cleaning
        pqPoint_i point_i;
        point_i.p = i;
        point_i.q = j;
        //std::cout << "not interpolated " << Grid[i][j] << std::endl;
        pqPoint point_interpolated = findMaximumInterpolated(point_i, Grid[i][j]);
        pqCollection.push_back(point_interpolated);
      }
    }
  }

}


pqPoint Retina::findMaximumInterpolated(pqPoint_i point_i, double w) {
  int p_i = point_i.p;
  int q_i = point_i.q;

  double p_mean = 0.;
  double q_mean = 0.;

  p_mean =  (pmin + pbinsize * (p_i - 0.5)) * Grid[p_i-1][q_i] +
            (pmin + pbinsize * (p_i + 0.5)) * Grid[p_i][q_i]   +
            (pmin + pbinsize * (p_i + 1.5)) * Grid[p_i+1][q_i];

  p_mean += (pmin + pbinsize * (p_i - 0.5)) * Grid[p_i-1][q_i-1] +
            (pmin + pbinsize * (p_i + 0.5)) * Grid[p_i][q_i-1]   +
            (pmin + pbinsize * (p_i + 1.5)) * Grid[p_i+1][q_i-1];

  p_mean += (pmin + pbinsize * (p_i - 0.5)) * Grid[p_i-1][q_i+1] +
            (pmin + pbinsize * (p_i + 0.5)) * Grid[p_i][q_i+1]   +
            (pmin + pbinsize * (p_i + 1.5)) * Grid[p_i+1][q_i+1];

  p_mean /= (Grid[p_i-1][q_i]   + Grid[p_i][q_i]   + Grid[p_i+1][q_i]   +
             Grid[p_i-1][q_i-1] + Grid[p_i][q_i-1] + Grid[p_i+1][q_i-1] +
             Grid[p_i-1][q_i+1] + Grid[p_i][q_i+1] + Grid[p_i+1][q_i+1]);


  q_mean =  (qmin + qbinsize * (q_i - 0.5)) * Grid[p_i][q_i-1] +
            (qmin + qbinsize * (q_i + 0.5)) * Grid[p_i][q_i]   +
            (qmin + qbinsize * (q_i + 1.5)) * Grid[p_i][q_i+1];

  q_mean += (qmin + qbinsize * (q_i - 0.5)) * Grid[p_i-1][q_i-1] +
            (qmin + qbinsize * (q_i + 0.5)) * Grid[p_i-1][q_i]   +
            (qmin + qbinsize * (q_i + 1.5)) * Grid[p_i-1][q_i+1];

  q_mean += (qmin + qbinsize * (q_i - 0.5)) * Grid[p_i+1][q_i-1] +
            (qmin + qbinsize * (q_i + 0.5)) * Grid[p_i+1][q_i]   +
            (qmin + qbinsize * (q_i + 1.5)) * Grid[p_i+1][q_i+1];

  q_mean /= (Grid[p_i][q_i-1]   + Grid[p_i][q_i]   + Grid[p_i][q_i+1]   +
             Grid[p_i-1][q_i-1] + Grid[p_i-1][q_i] + Grid[p_i-1][q_i+1] +
             Grid[p_i+1][q_i-1] + Grid[p_i+1][q_i] + Grid[p_i+1][q_i+1]);


  pqPoint point_o;
  point_o.p = p_mean;
  point_o.q = q_mean;
  point_o.w = w;
  //std::cout << "interpolated " << Grid[p_mean][j] << std::endl;
  return point_o;
}

void Retina::printMaxima() {
  unsigned int size = pqCollection.size();
  if ( view == XY )
    std::cout << "MAXIMA p , q , w (XY view):" << std::endl;
  else if ( view == RZ )
    std::cout << "MAXIMA p , q , w (RZ view):" << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    std::cout << pqCollection[i].p << " , " << pqCollection[i].q << " , " << pqCollection[i].w << std::endl;
  }
}

std::vector <pqPoint> Retina::getMaxima() {
  return pqCollection;
}

pqPoint Retina::getBestPQ() {
  unsigned int size = pqCollection.size();
  double max = 0.;
  pqPoint bestPQ;
  bestPQ.w = -1.;
  for (unsigned int i = 0; i < size; i++) {
    if (pqCollection[i].w > max) {
      max = pqCollection[i].w;
      bestPQ = pqCollection[i];
    }
  }
  //std::cout << " bestPQp " << bestPQ.p << "   bestPQq = " << bestPQ.q << " best w" << bestPQ.w << std::endl;
  return bestPQ;
}

