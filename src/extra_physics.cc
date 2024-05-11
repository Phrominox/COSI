// ----------------------------------------------------------------------------
// nexus | NexusPhysics.cc
//
// This class registers any new physics process defined in nexus.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "extra_physics.hh"

#include "wavelength_shifting.hh"

#include <G4OpticalPhoton.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessTable.hh>

  /// Macro that allows the use of this physics constructor
  /// with the generic physics list
  //G4_DECLARE_PHYSCONSTR_FACTORY(NexusPhysics);



  ExtraPhysics::ExtraPhysics():
    G4VPhysicsConstructor("ExtraPhysics")
  {
  }

  ExtraPhysics::~ExtraPhysics()
  {
  }

  void ExtraPhysics::ConstructParticle()
  {
    G4OpticalPhoton::Definition();
    //G4OpticalPhoton::OpticalPhotonDefinition();
  }

  void ExtraPhysics::ConstructProcess()
  {
    G4ProcessManager* pmanager = 0;

    // Add our own wavelength shifting process for the optical photon
    pmanager = G4OpticalPhoton::Definition()->GetProcessManager();
    if (!pmanager) {
      G4Exception("[ExtraPhysics]", "ConstructProcess()", FatalException,
        "G4OpticalPhoton without a process manager.");
    }
    WavelengthShifting* wls = new WavelengthShifting();
    pmanager->AddDiscreteProcess(wls);
  }