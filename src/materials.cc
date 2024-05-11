#include "materials.hh"

//#include "nain4.hh"
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-utils.hh"
//#include "n4-volumes.hh"

G4MaterialPropertiesTable* PTFE_props(double refl)
{
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
  // REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
//  std::vector<G4double> REFLECTIVITY = {0.8, 0.8};
  std::vector<G4double> REFLECTIVITY = {refl, refl};
  //std::vector<G4double> REFLECTIVITY = {1, 1};
  // REFLEXION BEHAVIOR
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  std::vector<G4double> specularlobe  = {0., 0.};
  // specular reflection about the average normal
  std::vector<G4double> specularspike = {0., 0.};
  // 180 degrees reflection.
  std::vector<G4double> backscatter   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  return n4::material_properties()
         .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
         .add("SPECULARLOBECONSTANT", ENERGIES, specularlobe)
         .add("SPECULARSPIKECONSTANT",ENERGIES, specularspike)
         .add("BACKSCATTERCONSTANT",  ENERGIES, backscatter)
         //.add("RINDEX", ENERGIES, {1.5, 1.5}) // if it doesn't reflect, it is just killed
         .done();
}

G4MaterialPropertiesTable* Vikuiti_props(double refl)
{
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
  // REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
//  std::vector<G4double> REFLECTIVITY = {0.8, 0.8};
  std::vector<G4double> REFLECTIVITY = {refl, refl};
  //std::vector<G4double> REFLECTIVITY = {1, 1};
  // REFLEXION BEHAVIOR
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  std::vector<G4double> specularlobe  = {0.5, 0.5};
  // specular reflection about the average normal
  std::vector<G4double> specularspike = {0.5, 0.5};
  // 180 degrees reflection.
  std::vector<G4double> backscatter   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  return n4::material_properties()
         .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
         .add("SPECULARLOBECONSTANT", ENERGIES, specularlobe)
         .add("SPECULARSPIKECONSTANT",ENERGIES, specularspike)
         .add("BACKSCATTERCONSTANT",  ENERGIES, backscatter)
         //.add("RINDEX", ENERGIES, {1.5, 1.5}) // if it doesn't reflect, it is just killed
         .done();
}

G4MaterialPropertiesTable* LAr_props()
{
    
    std::vector<G4double> energies;
    std::vector<G4double> wavelength;
    std::vector<G4double> rindex;
    std::vector<G4double> rayleigh;
    G4double wav;
    G4double rval;
    G4double lval;
    
    const double a_0 = 1.26;
    const double a_UV = 0.23;
    const double a_IR = 0.0023;
    const G4double gamma_UV = 106.6 * nm;
    const G4double gamma_IR = 908.3 * nm;
    
    for (double i=0.2;i<11.5;i=i+0.0113) { // 1000 steps
        energies.push_back(i*eV);
        wav = CLHEP::h_Planck * CLHEP::c_light / (i*eV);
        rval = sqrt(a_0 + (a_UV * pow(wav,2) / (pow(wav,2) - pow(gamma_UV,2))) + (a_IR * pow(wav,2) / (pow(wav,2) - pow(gamma_IR,2))));
        lval = (16*pow(CLHEP::pi,3)/(6*pow(wav,4))) * ((2.60943*pow(10,-24)*pow(1*cm,3))*pow((pow(rval,2) - 1)*(pow(rval,2) + 2)/3,2));
        
        wavelength.push_back(wav);
        rindex.push_back(rval); // refractive index
        rayleigh.push_back(pow(lval,-1)); // scattering length
    }
    
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    
//    for (size_t i = 100; i--;) 
//    G4cout << wavelength[i] / nm << "   " << rayleigh[i] / cm << G4endl;
    
    return n4::material_properties()
            .add("RINDEX", energies, rindex)
            .add("RAYLEIGH", energies, rayleigh)
            .add("ABSLENGTH", abs_energy, absLength)
            .done();
}

G4MaterialPropertiesTable* Al_props()
{
    
    std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> REFLECTIVITY = {0.7, 0.7};
    
    return n4::material_properties()
            .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
            .done();
}

G4MaterialPropertiesTable* Ti_props() // not a proper mapping of reflectivity, eventually digitize a plot
{
    std::vector<G4double> ENERGIES = {optPhotMinE_, 5*eV, 7.5*eV, optPhotMaxE_};
    std::vector<G4double> REFLECTIVITY = {0.6, 0.6, 0.8, 0.95};
    
    return n4::material_properties()
            .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
            .done();
}

G4MaterialPropertiesTable* tpb_props()
  {
    // Data from https://doi.org/10.1140/epjc/s10052-018-5807-z
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> TPB_rIndex      = {1.67    , 1.67};
    mpt->AddProperty("RINDEX", rIndex_energies, TPB_rIndex);

    // ABSORPTION LENGTH
    // Assuming no absorption except WLS
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH (Version NoSecWLS)
    // The NoSecWLS is forced by setting the WLS_absLength to noAbsLength_
    // for wavelengths higher than 380 nm where the WLS emission spectrum starts.
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      CLHEP::h_Planck * CLHEP::c_light / (380. * nm),  CLHEP::h_Planck * CLHEP::c_light / (370. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (360. * nm),  CLHEP::h_Planck * CLHEP::c_light / (330. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (320. * nm),  CLHEP::h_Planck * CLHEP::c_light / (310. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (300. * nm),  CLHEP::h_Planck * CLHEP::c_light / (270. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (250. * nm),  CLHEP::h_Planck * CLHEP::c_light / (230. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (210. * nm),  CLHEP::h_Planck * CLHEP::c_light / (190. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (170. * nm),  CLHEP::h_Planck * CLHEP::c_light / (150. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_,                 // ~6200 nm
      noAbsLength_,   50. * nm,     // 380 , 370 nm
      30. * nm,      30. * nm,      // 360 , 330 nm
      50. * nm,      80. * nm,      // 320 , 310 nm
      100. * nm,     100. * nm,     // 300 , 270 nm
      400. * nm,     400. * nm,     // 250 , 230 nm
      350. * nm,     250. * nm,     // 210 , 190 nm
      350. * nm,     400. * nm,     // 170 , 150 nm
      400. * nm                     // ~108 nm
    };

    //for (int i=0; i<WLS_abs_energy.size(); i++)
    //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    // Implemented with formula (7), with parameter values in table (3)
    // Sampling from ~380 nm to 600 nm <--> from 2.06 to 3.26 eV
    const G4int WLS_emi_entries = 120;
    std::vector<G4double> WLS_emi_energy;
    for (int i=0; i<WLS_emi_entries; i++)
      WLS_emi_energy.push_back(2.06 * eV + 0.01 * i * eV);

    std::vector<G4double> WLS_emiSpectrum;
    G4double A      = 0.782;
    G4double alpha  = 3.7e-2;
    G4double sigma1 = 15.43;
    G4double mu1    = 418.10;
    G4double sigma2 = 9.72;
    G4double mu2    = 411.2;

    for (int i=0; i<WLS_emi_entries; i++) {
      G4double wl = (CLHEP::h_Planck * CLHEP::c_light / WLS_emi_energy[i]) / nm;
      WLS_emiSpectrum.push_back(A * (alpha/2.) * exp((alpha/2.) *
                          (2*mu1 + alpha*pow(sigma1,2) - 2*wl)) *
                          erfc((mu1 + alpha*pow(sigma1,2) - wl) / (sqrt(2)*sigma1)) +
                          (1-A) * (1 / sqrt(2*pow(sigma2,2)*3.1416)) *
                                exp((-pow(wl-mu2,2)) / (2*pow(sigma2,2))));
      // G4cout << "* TPB WLSemi:  " << std::setw(4)
      //        << wl << " nm -> " << WLS_emiSpectrum[i] << G4endl;
    };
    mpt->AddProperty("WLSCOMPONENT", WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1.2 * ns);

    // WLS Quantum Efficiency
    // According to the paper, the QE of TPB depends on the incident wavelength.
    // As Geant4 doesn't allow this possibility, it is set to the value corresponding
    // to Xe scintillation spectrum peak.
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.65);

    return mpt;
  }

G4Material* TPB()
{
  G4String name = "TPB"; // Tetraphenyl butadiene
  G4Material* mat = G4Material::GetMaterial(name, false);
  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();
    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    mat = new G4Material(name, 1*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 22);
    mat->AddElement(C, 28);
  }
  return mat;
}

G4Material* tpb_with_props(){
  auto tpb = TPB();
  tpb -> SetMaterialPropertiesTable(tpb_props());
  return tpb;
}