// code draws random samples from the LAr volume, generates X scintillation photons at each point, then saves out all necessary info like energy, final position on an SiPM (for later being able to do cuts on filling factor)

// so assume each unique event ID has X trials, can calculate an efficiency of reaching a detector for each group of wavelength prominent in the problem (for this one 128 nm and 400-500 nm) for which the SiPM has different QE's
// total efficiency is the sum, though should be weighted by QE

#include <n4-all.hh>
#include <n4-boolean-shape.hh>

#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4IntersectionSolid.hh>

#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include "G4AnalysisManager.hh"
#include <G4ThreeVector.hh>

#include <G4GenericMessenger.hh>

#include <G4PrimaryParticle.hh>
#include <G4String.hh>
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include <G4Event.hh>           // needed to inject primary particles into an event
#include <G4Box.hh>             // for creating shapes in the geometry
#include <G4Sphere.hh>          // for creating shapes in the geometry
#include <FTFP_BERT.hh>         // our choice of physics list
#include <G4RandomDirection.hh> // for launching particles in random directions

#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4IonTable.hh>

#include "extra_physics.hh"
#include "materials.hh"

#include <cstdlib>
#include <fstream>
using std::ofstream;
#include <iostream>

ofstream output;
G4ThreeVector event_initial;

n4::sensitive_detector* sensor(unsigned& n_event) {
  auto process_hits = [&] (G4Step* step) {
    
    step -> GetTrack() -> SetTrackStatus(fStopAndKill);
    
    auto pre           = step -> GetPreStepPoint();
//    auto post           = step -> GetPostStepPoint();
    auto photon_energy = pre  -> GetKineticEnergy();
    
    auto pos = pre -> GetPosition();
    auto pos_initial = step -> GetTrack() -> GetVertexPosition(); // gets where it was converted, if done so
    
    auto cell = pre -> GetTouchable() -> GetVolume() -> GetName();
   
//    G4cout << "E: " << photon_energy << "     pos: " << pos << G4endl;
    
    output << n_event << " " << event_initial << " " << cell << " " << photon_energy << " " << pos << "\n";
    
    return true; // See https://jacg.github.io/nain4/explanation/process-hits-return-value.html
  };
  return new n4::sensitive_detector{"detector", process_hits};
}


auto geometry(unsigned& n_event) {
    
    double ptfe_thick = 2*mm;
    double tpb_thick = 1*um;
    double cryo_rad = 19.9*cm; // inner radius of the cryostat
    double cryo_len = 50*cm; // inner total length of the cryostat containing the veto + detectors
    double wall_thick = 3*mm; // cryostat wall thickness
    double vacuum_buffer = 3*cm;
    
    auto air = n4::material("G4_AIR");
    auto vac = n4::material("G4_AIR");
    auto copper = n4::material("G4_Cu");
    auto Al = n4::material("G4_Al");
    auto steel  = n4::material("G4_STAINLESS-STEEL");
    auto ptfe   = n4::material("G4_TEFLON");
    auto csi = n4::material("G4_CESIUM_IODIDE");
    auto LAr = n4::material("G4_lAr");
    auto Ti = n4::material("G4_Ti");
    auto tpb = tpb_with_props();
    
    air->SetMaterialPropertiesTable(n4::material_properties().add("RINDEX", {optPhotMinE_, optPhotMaxE_}, {1.0003, 1.0003}).done());
    LAr->SetMaterialPropertiesTable(LAr_props());
    // copper->SetMaterialPropertiesTable(n4::material_properties().add("RINDEX", {optPhotMinE_, optPhotMaxE_}, {1.0003, 1.0003}).done());
    
 //   Al -> SetMaterialPropertiesTable(Al_props()); // obviated by making a logical surface?
 //   Ti -> SetMaterialPropertiesTable(Ti_props());
    
    double refl_ptfe = 0.98;//0.98; // also true of Vikuiti @ 98%
    G4OpticalSurface* ptfe_reflector = new G4OpticalSurface("ptfe_reflector", unified, polished, dielectric_metal);
    ptfe_reflector->SetMaterialPropertiesTable(PTFE_props(refl_ptfe));
    
    // don't end up using this anymore, since holder is now Cu and covered by waveshifter
    G4OpticalSurface* Al_reflector = new G4OpticalSurface("Al_reflector", unified, polished, dielectric_metal);
    Al_reflector->SetMaterialPropertiesTable(Al_props());
    
    G4OpticalSurface* Ti_reflector = new G4OpticalSurface("Ti_reflector", unified, polished, dielectric_metal);
    Ti_reflector->SetMaterialPropertiesTable(Ti_props());
    
    auto world = n4::box("world").xyz(1*m,1*m,1*m).volume(air);

    auto outer_cryostat = n4::tubs("outer_cryostat").r(cryo_rad + wall_thick + vacuum_buffer + wall_thick).z(cryo_len + wall_thick*2 + vacuum_buffer*2 + wall_thick*2).place(Ti).in(world).now();
    auto vacuum = n4::tubs("cryostat").r(cryo_rad + wall_thick + vacuum_buffer).z(cryo_len + wall_thick*2 + vacuum_buffer*2).place(vac).in(outer_cryostat).now();
    
    auto cryostat = n4::tubs("cryostat").r(cryo_rad + wall_thick).z(cryo_len + wall_thick*2).place(Ti).in(vacuum).now();
    
    auto reflector_cyl = n4::tubs("reflector_cyl").r((cryo_rad)).r_delta(ptfe_thick).z(cryo_len).place(ptfe).in(cryostat).now();
    auto reflector_cap_1 = n4::tubs("reflector_cap_1").r((cryo_rad)).z(ptfe_thick).sub(n4::box("detector_hole_refl1").xyz(20*cm,20*cm,ptfe_thick)).place(ptfe).in(cryostat).at_z((cryo_len/2 - ptfe_thick/2)).now();
    auto reflector_cap_2 = n4::tubs("reflector_cap_2").r((cryo_rad)).z(ptfe_thick).sub(n4::box("detector_hole_refl2").xyz(20*cm,20*cm,ptfe_thick)).place(ptfe).in(cryostat).at_z(-(cryo_len/2 - ptfe_thick/2)).now();
    
    auto tpb_vol = n4::tubs("tpb").r(cryo_rad - ptfe_thick).r_delta(tpb_thick).z(cryo_len).volume(tpb);
    auto tpb_coat = n4::place(tpb_vol).in(cryostat).now();

    auto tpb_cap_vol_1 = n4::tubs("tpb_cap_1").r((cryo_rad)).z(tpb_thick).sub(n4::box("detector_hole_cap1").xyz(20*cm,20*cm,tpb_thick)).volume(tpb);
    auto tpb_cap_coat_1 = n4::place(tpb_cap_vol_1).in(cryostat).at_z((cryo_len/2 - ptfe_thick - tpb_thick/2)).now();
    auto tpb_cap_vol_2 = n4::tubs("tpb_cap_2").r((cryo_rad)).z(tpb_thick).sub(n4::box("detector_hole_cap2").xyz(20*cm,20*cm,tpb_thick)).volume(tpb);
    auto tpb_cap_coat_2 = n4::place(tpb_cap_vol_2).in(cryostat).at_z(-(cryo_len/2 - ptfe_thick - tpb_thick/2)).now();
    
    auto liquid = n4::tubs("liquid").r(cryo_rad - ptfe_thick - tpb_thick).z(cryo_len - ptfe_thick*2 - tpb_thick*2).place(LAr).in(cryostat).now();
    
    // new G4LogicalBorderSurface("Ti_reflector", liquid, cryostat, Ti_reflector); // meant to be the endcap of the cryostat, can also cover this with waveshifter

    auto reflector_plane_1a = n4::box("reflector_plane_1a").xyz(221*mm,154*mm,ptfe_thick);
    auto reflector_plane_1b = n4::box("reflector_plane_1b").xyz(77*mm,239*mm,ptfe_thick);
    auto reflector_plane_1 = reflector_plane_1a.add(reflector_plane_1b).place(ptfe).in(liquid).at_z((16*cm + 1*cm + ptfe_thick/2)).now();
    auto reflector_plane_2a = n4::box("reflector_plane_2a").xyz(221*mm,154*mm,ptfe_thick);
    auto reflector_plane_2b = n4::box("reflector_plane_2b").xyz(77*mm,239*mm,ptfe_thick);
    auto reflector_plane_2 = reflector_plane_2a.add(reflector_plane_2b).place(ptfe).in(liquid).at_z(-(16*cm + 1*cm + ptfe_thick/2)).now();

    auto tpb_plane_vol_1a = n4::box("tpb_plane_1a").xyz(221*mm,154*mm,tpb_thick); // means photons between plane and holder/crystals go from LAr to refl and are killed anyway
    auto tpb_plane_vol_1b = n4::box("tpb_plane_1b").xyz(77*mm,239*mm,tpb_thick);
    auto tpb_plane_vol_1 = tpb_plane_vol_1a.add(tpb_plane_vol_1b).volume(tpb);
    auto tpb_plane_coat_1 = n4::place(tpb_plane_vol_1).in(liquid).at_z((16*cm + 1*cm + ptfe_thick + tpb_thick/2)).now();
    auto tpb_plane_vol_2a = n4::box("tpb_plane_2a").xyz(221*mm,154*mm,tpb_thick);
    auto tpb_plane_vol_2b = n4::box("tpb_plane_2b").xyz(77*mm,239*mm,tpb_thick);
    auto tpb_plane_vol_2 = tpb_plane_vol_2a.add(tpb_plane_vol_2b).volume(tpb);
    auto tpb_plane_coat_2 = n4::place(tpb_plane_vol_2).in(liquid).at_z(-(16*cm + 1*cm + ptfe_thick + tpb_thick/2)).now();
    
    G4OpticalSurface* tpb_surf_cyl = new G4OpticalSurface("tpb_opsurf_cyl", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("tpb_opsurf_cyl", tpb_vol, tpb_surf_cyl);

    G4OpticalSurface* tpb_surf_cap_1 = new G4OpticalSurface("tpb_opsurf_cap_1", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("tpb_opsurf_cap_1", tpb_cap_vol_1, tpb_surf_cap_1);
    G4OpticalSurface* tpb_surf_cap_2 = new G4OpticalSurface("tpb_opsurf_cap_2", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("tpb_opsurf_cap_2", tpb_cap_vol_2, tpb_surf_cap_2);

    G4OpticalSurface* tpb_surf_plane_1 = new G4OpticalSurface("tpb_opsurf_plane_1", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("tpb_opsurf_plane_1", tpb_plane_vol_1, tpb_surf_plane_1);
    G4OpticalSurface* tpb_surf_plane_2 = new G4OpticalSurface("tpb_opsurf_plane_2", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("tpb_opsurf_plane_2", tpb_plane_vol_2, tpb_surf_plane_2);
    
    // reflecting emmisions that pass through the TPB
    new G4LogicalBorderSurface("internal_reflector_1", tpb_coat, reflector_cyl, ptfe_reflector);
    new G4LogicalBorderSurface("internal_reflector_2", tpb_cap_coat_1, reflector_cap_1, ptfe_reflector);
    new G4LogicalBorderSurface("internal_reflector_3", tpb_cap_coat_2, reflector_cap_2, ptfe_reflector);
    new G4LogicalBorderSurface("internal_reflector_4", tpb_plane_coat_1, reflector_plane_1, ptfe_reflector);
    new G4LogicalBorderSurface("internal_reflector_5", tpb_plane_coat_2, reflector_plane_2, ptfe_reflector);

//  material is currently air because index of refraction is already set and the particle is killed anyway
    auto sipm1 = n4::box("sipm1").xyz(20*cm,20*cm,5*mm).sensitive(sensor(n_event)).place(air).in(liquid).at_z(cryo_len/2 - 5*mm/2).now();
    auto sipm2 = n4::box("sipm2").xyz(20*cm,20*cm,5*mm).sensitive(sensor(n_event)).place(air).in(liquid).at_z(-(cryo_len/2 - 5*mm/2)).now();
    // it is not necessary to use a sensitive detector and one can just use boolean conditions in the "step" action to define when one reaches the correct volume, but for now it doesn't matter
    
    auto tpb_sup_vol_1 = n4::box("tpb_sup_1").xyz(221*mm + ptfe_thick + tpb_thick,154*mm + ptfe_thick + tpb_thick,32*cm - 1*mm).sub(n4::box("hole1").xyz(221*mm + ptfe_thick,154*mm + ptfe_thick,32*cm)).sub(n4::box("hole2").xyz(77*mm + ptfe_thick + tpb_thick,239*mm + ptfe_thick + tpb_thick,32*cm));
    auto tpb_sup_vol_2 = n4::box("tpb_sup_2").xyz(77*mm + ptfe_thick + tpb_thick,239*mm + ptfe_thick + tpb_thick,32*cm - 1*mm).sub(n4::box("hole3").xyz(77*mm + ptfe_thick,239*mm + ptfe_thick,32*cm)).sub(n4::box("hole4").xyz(221*mm + ptfe_thick + tpb_thick,154*mm + ptfe_thick + tpb_thick,32*cm));
    auto tpb_sup_vol = tpb_sup_vol_1.add(tpb_sup_vol_2).volume(tpb);
    auto tpb_sup_coat = n4::place(tpb_sup_vol).in(liquid).now();

    auto reflector_sup_1 = n4::box("refl_sup_1").xyz(221*mm + ptfe_thick,154*mm + ptfe_thick,32*cm - 1*mm).sub(n4::box("hole5").xyz(221*mm,154*mm,32*cm)).sub(n4::box("hole6").xyz(77*mm + ptfe_thick,239*mm + ptfe_thick,32*cm));
    auto reflector_sup_2 = n4::box("refl_sup_2").xyz(77*mm + ptfe_thick,239*mm + ptfe_thick,32*cm - 1*mm).sub(n4::box("hole7").xyz(77*mm,239*mm,32*cm)).sub(n4::box("hole8").xyz(221*mm + ptfe_thick,154*mm + ptfe_thick,32*cm));
    auto reflector_sup = reflector_sup_1.add(reflector_sup_2).place(ptfe).in(liquid).now();

    auto holder_1 = n4::box("holder").xyz(221*mm,154*mm,1*cm);
    auto holder_2 = n4::box("holder2").xyz(77*mm,239*mm,1*cm);    

    // reduced reflection of LAr emission, but decent for TPB emission
    //new G4LogicalBorderSurface("Al_reflector", liquid, holder, Al_reflector);

    G4OpticalSurface* tpb_surf_sup = new G4OpticalSurface("tpb_opsurf_sup", glisur, ground,
                           dielectric_dielectric, .01);
    new G4LogicalSkinSurface("tpb_opsurf_sup", tpb_sup_vol, tpb_surf_sup);

    new G4LogicalBorderSurface("internal_reflector_6", tpb_sup_coat, reflector_sup, ptfe_reflector);
    
    // crystals just stop photons at this time (without refractive index defined), will be covered by APD's
    auto crystal_box = n4::box("crystal").xyz(6.7*cm,6.7*cm,32*cm);
    auto crystal_volume = n4::volume(crystal_box.solid(), csi);
    auto crystal = n4::place(crystal_volume).in(liquid);
    crystal.clone().copy_no(1).now();
    crystal.clone().copy_no(2).at_y(81*mm).now();
    crystal.clone().copy_no(3).at_x(72*mm).at_y(38.5*mm).now();
    crystal.clone().copy_no(4).at_x(72*mm).at_y(-38.5*mm).now();
    crystal.clone().copy_no(5).at_y(-81*mm).now();
    crystal.clone().copy_no(6).at_x(-72*mm).at_y(-38.5*mm).now();
    crystal.clone().copy_no(7).at_x(-72*mm).at_y(38.5*mm).now();

    auto no_tran = G4ThreeVector(0*mm,0*mm,0*mm);
    G4RotationMatrix* no_rot = new G4RotationMatrix();

    auto t2 = G4ThreeVector(0*mm,81*mm,0*mm); // same as the translations of the individual crystals above
    auto t3 = G4ThreeVector(72*mm,38.5*mm,0*mm);
    auto t4 = G4ThreeVector(72*mm,-38.5*mm,0*mm);
    auto t5 = G4ThreeVector(0*mm,-81*mm,0*mm);
    auto t6 = G4ThreeVector(-72*mm,-38.5*mm,0*mm);
    auto t7 = G4ThreeVector(-72*mm,38.5*mm,0*mm);

    G4String label = "composite_holder";
    auto composite_holder = new G4UnionSolid{label,holder_1.solid(),holder_2.solid(),no_rot,no_tran};

    G4String label_crys1 = "holder_minus_crys1";
    auto holder_minus_crys1 = new G4SubtractionSolid{label_crys1,composite_holder,crystal_box.solid(),no_rot,no_tran};
    G4String label_crys2 = "holder_minus_crys2";
    auto holder_minus_crys2 = new G4SubtractionSolid{label_crys2,holder_minus_crys1,crystal_box.solid(),no_rot,t2};
    G4String label_crys3 = "holder_minus_crys3";
    auto holder_minus_crys3 = new G4SubtractionSolid{label_crys3,holder_minus_crys2,crystal_box.solid(),no_rot,t3};
    G4String label_crys4 = "holder_minus_crys4";
    auto holder_minus_crys4 = new G4SubtractionSolid{label_crys4,holder_minus_crys3,crystal_box.solid(),no_rot,t4};
    G4String label_crys5 = "holder_minus_crys5";
    auto holder_minus_crys5 = new G4SubtractionSolid{label_crys5,holder_minus_crys4,crystal_box.solid(),no_rot,t5};
    G4String label_crys6 = "holder_minus_crys6";
    auto holder_minus_crys6 = new G4SubtractionSolid{label_crys6,holder_minus_crys5,crystal_box.solid(),no_rot,t6};
    G4String label_crys7 = "holder_minus_crys7";
    auto holder_minus_crys7 = new G4SubtractionSolid{label_crys7,holder_minus_crys6,crystal_box.solid(),no_rot,t7};

    auto holder_volume = n4::volume(holder_minus_crys7,copper);
    auto holder = n4::place(holder_volume).in(liquid);

    holder.clone().copy_no(1).at_z(-15.5*cm).now();
    holder.clone().copy_no(2)               .now();
    holder.clone().copy_no(3).at_z(15.5*cm).now();

    return n4::place(world).now();
}


auto generator() {
return [&](G4Event* event) {
    auto particle_type = n4::find_particle("opticalphoton");
//    auto particle = new G4PrimaryParticle(particle_type);
    
    
    auto energy = 9.68627 * eV; // LAr emission at 128 nm
    
    do {
        
        double r = 20.2*cm * sqrt(G4UniformRand());
        double z = G4UniformRand() * 50*cm - 25*cm; // to generate throughout volume if center is origin
        double phi = G4UniformRand() * 2*M_PI;
        
        double x = r * cos(phi);
        double y = r * sin(phi);
        
        // based on crystal holder dimensions with 5 mm of Cu between liquid and CsI
        if ((std::abs(x) < 110.5*mm && std::abs(y) < 77*mm) || (std::abs(x) < 38.5*mm && std::abs(y) < 119.5*mm)) {
            if (std::abs(z) < (32*cm/2 + 1*cm + 2*mm + 1*um)) { // need to make sure the sample outside of the additonal waveshifter planes on either end (=> + 1 cm + ptfe_thick + tpb_thick)
                continue;
            }
        }
        
        auto initial_pos = G4ThreeVector(x, y, z);
    
        int sample_size = 1; // how many photons should be emitted from each sample point in the volume, then i don't need to bin them volumetrically and can smooth over the nework of calculated efficiencies to get a heat map
    
        for (int i=0; i<sample_size; i++) {
            
            auto particle = new G4PrimaryParticle(particle_type);
            
            auto p = G4RandomDirection();
            auto po = G4RandomDirection();
            
            particle -> SetMomentumDirection(p);
            particle -> SetPolarization(po);
            particle -> SetKineticEnergy(energy);
            
            auto vertex   = new G4PrimaryVertex(initial_pos, 0);
            vertex -> SetPrimary(particle);
            event  -> AddPrimaryVertex(vertex);
        }
        
        break;
        
    } while (true);
    
    };
}

n4::actions* create_actions(unsigned& n_event) {
    
    auto per_step = [&] (const G4Step* step) {
//        auto cell = step -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName();
//        if (cell == "detector") {
//           auto pos = step -> GetPreStepPoint() -> GetPosition();
//            std::cout << cell << " " << pos << std::endl;
        
//        auto scoring_vol = n4::find_logical("detector");
//        auto current_vol = step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume();

        //auto scoring_vol = n4::find_logical("sipm");
        auto current_vol = step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetName();

        
//        std::cout << scoring_vol << std::endl;
//        std::cout << current_vol << std::endl;
//        std::cout << next_vol << std::endl;
        
        if (current_vol == "sipm") {
            auto pos = step -> GetPreStepPoint() -> GetPosition();
            //std::cout << scoring_vol << " " << pos << std::endl;
        }
    };
    
    auto per_event_beg = [&] (const G4Event* event) {
        event_initial = event -> GetPrimaryVertex() -> GetPosition();
    };
    
    auto per_event_end = [&] (const G4Event*) {
       std::cout << "end of event " << n_event << std::endl;
       n_event++;
    };
    
    return (new n4::        actions{generator()  })
 -> set( (new n4::   event_action{                  })
        -> begin(per_event_beg)
        -> end(per_event_end) )
 -> set(  new n4::stepping_action{per_step});
}

int main(int argc, char** argv)
{
    unsigned n_event = 0;
    
    G4int verbosity=0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list -> ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    physics_list -> RegisterPhysics(new G4RadioactiveDecayPhysics);
    physics_list -> RegisterPhysics(new ExtraPhysics{});
    
    output.open ("caca.txt");
    
    n4::run_manager::create()	
    .ui("LAr", argc, argv)
    .macro_path("macs")
 //   .apply_cli_early()
    .physics(physics_list)
    .geometry([&] {return geometry(n_event);})
    .actions(create_actions(n_event))
 //   .apply_cli_late()
    .run();
    
    output.close();
}

