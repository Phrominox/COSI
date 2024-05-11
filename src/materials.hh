
#include "G4MaterialPropertiesTable.hh"
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include "G4Material.hh"
#include "cmath"
#include <CLHEP/Units/PhysicalConstants.h>

#include <n4-all.hh>

const G4double optPhotMinE_   =  0.2  * eV;
const G4double optPhotMaxE_   = 11.5  * eV;
const G4double noAbsLength_   = 1.e8  * m;

G4MaterialPropertiesTable* PTFE_props(double refl);
G4MaterialPropertiesTable* Vikuiti_props(double refl);
G4MaterialPropertiesTable* LAr_props();
G4MaterialPropertiesTable* Al_props();
G4MaterialPropertiesTable* Ti_props();
G4MaterialPropertiesTable* tpb_props();

G4Material* TPB();
G4Material* tpb_with_props();
                                          
