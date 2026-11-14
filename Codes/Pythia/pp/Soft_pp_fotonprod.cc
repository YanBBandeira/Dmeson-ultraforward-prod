#include <iostream>
#include "Pythia8/Pythia.h"
#include <string>
#include <filesystem>

using namespace Pythia8;
namespace fs = std::filesystem;

// Função principal
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Uso: " << argv[0] << "<ecm_em_GeV>" << std::endl;
  }

  std::string eCM_str = "13000";
  Pythia pythia;

  pythia.readString("Beams:idA = 2212"); // proton
  pythia.readString("Beams:idB = 2212"); // proton
  pythia.readString("Beams:eCM = " + eCM_str);

  //pythia.readString("SoftQCD:all = on");
  pythia.readString("SoftQCD:nonDiffractive = on");
  //pythia.readString("MultipartonInteractions:pT0Ref = 3"); 

  //FLAGS PARA EPA
  pythia.readString("Photon:ProcessType = 0");
  pythia.readString("PhotonParton:all = on");

 
  pythia.readString("PhaseSpace:pThatMin = 0."); // Minimal pT cut
  pythia.readString("Photon:Wmin = 5.");
  pythia.readString("Photon:sampleQ2 = off");
  pythia.readString("PDF:beamA2gamma = on ");
  pythia.readString("PDF:proton2gammaSet = 2"); // Use a virtuality-dependent flux similar to the approximation by Drees and Zeppenfeld. Uses a dipole form factor to appropriately suppress high virtualities.
             
  pythia.readString("411:mayDecay = off");  // D+
  pythia.readString("-411:mayDecay = off"); // D-
  pythia.readString("421:mayDecay = off");  // D0
  pythia.readString("-421:mayDecay = off"); // D0bar
  pythia.readString("413:mayDecay = off");  // D*+
  pythia.readString("-413:mayDecay = off"); // D*-
  pythia.readString("423:mayDecay = off");  // D*0
  pythia.readString("-423:mayDecay = off"); // D*0bar
  pythia.readString("415:mayDecay = off");  // D2*+
  pythia.readString("-415:mayDecay = off"); // D2*-
  pythia.readString("425:mayDecay = off");  // D2*0
  pythia.readString("-425:mayDecay = off"); // D2*0bar
  pythia.readString("431:mayDecay = off");  // Ds+
  pythia.readString("-431:mayDecay = off"); // Ds-
  pythia.readString("433:mayDecay = off");  // Ds*+
  pythia.readString("-433:mayDecay = off"); // Ds*-
  pythia.readString("435:mayDecay = off");  // Ds2*
  pythia.readString("-435:mayDecay = off"); // Ds2*bar

  if (!pythia.init())
  {
    cout << "Erro: falha na inicialização do Pythia." << endl;
    return 1;
  }
  cout << "Pythia inicializado com sucesso." << endl;


  double eta_min = 0.0;
  double eta_max = 10.0;
  double y_min = 0.0; 
  double y_max = 10.0;
  double pt_min = 0.0;
  double pt_max =  50.0;

  int bins = 100;
  double binWidth_eta = (eta_max - eta_min) / bins;
  double binWidth_y   = (y_max - y_min)/ bins;
  double binWidth_pt  = (pt_max - pt_min) / bins;

  Hist dNeta_D0meson_eta = Hist(" ", bins, eta_min, eta_max);
  Hist dNeta_Dsmeson_eta = Hist(" ", bins, eta_min, eta_max);

  Hist dNeta_D0meson_y = Hist(" ", bins, y_min, y_max);
  Hist dNeta_Dsmeson_y = Hist(" ", bins, y_min, y_max);

  Hist dNeta_D0meson_pt = Hist(" ", bins, pt_min, pt_max);
  Hist dNeta_Dsmeson_pt = Hist(" ", bins, pt_min, pt_max);


  Hist mult_D0meson_eta = Hist(" ", 100, eta_min,eta_max, 99.5);
  Hist mult_Dsmeson_eta = Hist(" ", 100, eta_min,eta_max, 99.5);



  int nEvent = 10000000; // 1000000;
  int nAcceptedEvents = 0;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythia.next())
      continue;

    nAcceptedEvents++;
    int countD0 = 0;
    int countDs = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    {

      int id = pythia.event[i].id();

      if (pythia.event[i].isFinal() && abs(pythia.event[i].id()) == 421) //D0 meson
      {
        double eta = pythia.event[i].eta();
        double pt = pythia.event[i].pT();
        double y = pythia.event[i].y();
        if (eta >= eta_min && eta <= eta_max)
        {
          dNeta_D0meson_eta.fill(eta);
          countD0++;
        }
        if (y >= y_min && y <= y_max){
          dNeta_D0meson_y.fill(y);
        }
        if (pt >= pt_min && pt <= pt_max){
          dNeta_D0meson_pt.fill(pt);
        }
      }

      if (pythia.event[i].isFinal() && abs(pythia.event[i].id()) == 431) //Ds meson
      {
        double eta = pythia.event[i].eta();
        double pt = pythia.event[i].pT();
        double y = pythia.event[i].y();
        if (eta >= eta_min && eta <= eta_max)
        {
          dNeta_Dsmeson_eta.fill(eta);
          countDs++;
        }
        if (y >= y_min && y <= y_max){
          dNeta_Dsmeson_y.fill(y);
        }
        if (pt >= pt_min && pt <= pt_max){
          dNeta_Dsmeson_pt.fill(pt);
        }
      }
    }
    mult_D0meson_eta.fill(countD0);
    mult_Dsmeson_eta.fill(countDs);
  }

  // NOrmalização
  double normFactor_eta = 1.0 / (nAcceptedEvents * binWidth_eta);
  double normFactor_y   = 1.0 / (nAcceptedEvents * binWidth_y);
  double normFactor_pt  = 1.0 / (nAcceptedEvents * binWidth_pt);


   // 2. Criação das distribuições de Seção de Choque (dsigma/dx)
  // pythia.info.sigmaGen() retorna a seção de choque total em milibarns (mb)
  double sigmaGen_mb = pythia.info.sigmaGen(); 
  
  // Opcional: multiplique por 1e3 para microbarns (ub) ou 1e9 para picobarns (pb)
  double unitConv = 1.0; // mantém em mb

  Hist dSigma_D0meson_eta = dNeta_D0meson_eta * (sigmaGen_mb * unitConv);
  Hist dSigma_Dsmeson_eta = dNeta_Dsmeson_eta * (sigmaGen_mb * unitConv);

  Hist dSigma_D0meson_y   = dNeta_D0meson_y   * (sigmaGen_mb * unitConv);
  Hist dSigma_Dsmeson_y   = dNeta_Dsmeson_y   * (sigmaGen_mb * unitConv);

  Hist dSigma_D0meson_pt  = dNeta_D0meson_pt  * (sigmaGen_mb * unitConv);
  Hist dSigma_Dsmeson_pt  = dNeta_Dsmeson_pt  * (sigmaGen_mb * unitConv);

  // 3. Salvamento estruturado
  std::string baseOut = "../../Output/pp/Soft_pp/eCM_" + eCM_str + "GeV/";
  fs::create_directories(baseOut);

  // Salvando distribuições de N (1/Nev * dN/dx)
  dNeta_D0meson_eta.table(baseOut + "dNdeta_D0meson_eta.dat");
  dNeta_D0meson_y.table(baseOut + "dNdy_D0meson_y.dat");
  dNeta_D0meson_pt.table(baseOut + "dNdpt_D0meson_pt.dat");
  dNeta_Dsmeson_eta.table(baseOut + "dNdeta_Dsmeson_eta.dat");
  dNeta_Dsmeson_y.table(baseOut + "dNdy_Dsmeson_y.dat");
  dNeta_Dsmeson_pt.table(baseOut + "dNdpt_Dsmeson_pt.dat");

  // Salvando distribuições de Seção de Choque (dsigma/dx em mb)
  dSigma_D0meson_eta.table(baseOut + "dSigmadeta_D0meson_eta.dat");
  dSigma_D0meson_y.table(baseOut + "dSigmady_D0meson_y.dat");
  dSigma_D0meson_pt.table(baseOut + "dSigmadpt_D0meson_pt.dat");
  dSigma_Dsmeson_eta.table(baseOut + "dSigmadeta_Dsmeson_eta.dat");
  dSigma_Dsmeson_y.table(baseOut + "dSigmady_Dsmeson_y.dat");
  dSigma_Dsmeson_pt.table(baseOut + "dSigmadpt_Dsmeson_pt.dat");

  mult_D0meson_eta.table(baseOut + "mult_D0meson_eta.dat");
  mult_Dsmeson_eta.table(baseOut + "mult_Dsmeson_eta.dat");

  cout << "Secao de choque total gerada: " << sigmaGen_mb << " mb" << endl;
  cout << "Simulacao finalizada e arquivos salvos em " << baseOut << endl;
}
