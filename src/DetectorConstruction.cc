#include "DetectorConstruction.hh"

#include "EventStore.hh"
#include "CaloSD.hh"

// GeoModel builder
#include "MaterialManager.h"
#include "ConfigReader.h"
#include "CalorimeterBuilder.h"

#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/Units.h"

#include "GeoModelKernel/GeoNameTag.h"
#include "GeoModelKernel/GeoTransform.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>



// GeoModel2G4
#include "GeoModel2G4/Geo2G4AssemblyFactory.h"
#include "GeoModel2G4/Geo2G4AssemblyVolume.h"
#include "GeoModel2G4/VolumeBuilder.h"

// Geant4
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// Export to gdml
#include "G4GDMLParser.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include <cstdlib>


using namespace GeoModelKernelUnits;

static G4VisAttributes* MakeVis(double r, double g, double b, double a=1.0) {
    auto* v = new G4VisAttributes(G4Colour(r, g, b, a));
    v->SetVisibility(true);
    v->SetForceSolid(true);
    return v;
}

double GetSystemThickness(std::vector<int> vec_layers,double thickness_wide, double thickness_thin, double thickness_hpl, double thickness_passive, double thickness_airgap){

    double total_thickness = 0;

    for(int N_layer=0;N_layer<vec_layers.size();N_layer++){
        int layercode = vec_layers[N_layer];
        switch(layercode){
            case 1: total_thickness += thickness_wide; break;
            case 2: total_thickness += thickness_wide; break;
            case 3: total_thickness += thickness_thin; break;
            case 4: total_thickness += thickness_thin; break;
            case 5: total_thickness += thickness_hpl; break;
            case 6: total_thickness += thickness_hpl; break;
            case 7: total_thickness += thickness_passive; break;
            case 8: total_thickness += thickness_airgap; break;
        }
    }

    return total_thickness;
}


DetectorConstruction::DetectorConstruction(EventStore* store, std::string cfgFile, int w_g)
: m_store(store), m_cfgFile(std::move(cfgFile)) {write_gdml = w_g;}

DetectorConstruction::~DetectorConstruction() {}

GeoPhysVol* DetectorConstruction::buildGeoModelWorld()
{
  MaterialManager MM;
  auto* air = MM.air();

  auto cfg = readConfigFile(m_cfgFile);
  
  double aworldZ = GetSystemThickness(cfg.layers,cfg.scint_thickness_mm,cfg.scint_thickness_mm,cfg.hpl_thickness_mm,cfg.lead_thickness_mm,cfg.airgap_mm);
  aworldZ += cfg.gap_ecal_hcal_mm;
  aworldZ += GetSystemThickness(cfg.layers2,cfg.scint_thickness_mm,cfg.scint_thickness_mm,cfg.hpl_thickness_mm,cfg.iron_thickness_mm,cfg.airgap_mm);
  
    
  const int nx = std::max(1, cfg.module_nx);
  const int ny = std::max(1, cfg.module_ny);
  const double plateSize = cfg.plate_xy_mm * GeoModelKernelUnits::mm;
  const double airThickness = cfg.boundary_air_mm * GeoModelKernelUnits::mm;
  const double steelThickness = cfg.boundary_steel_mm * GeoModelKernelUnits::mm;
  if (plateSize <= 0.0 || airThickness <= 0.0 || steelThickness <= 0.0)
    throw std::runtime_error("Module and casing dimensions must be positive");

  const double casingSize = plateSize + 2.0 * (airThickness + steelThickness);
  const double pitchX = (cfg.module_pitch_x_mm > 0.0
      ? cfg.module_pitch_x_mm * GeoModelKernelUnits::mm : casingSize);
  const double pitchY = (cfg.module_pitch_y_mm > 0.0
      ? cfg.module_pitch_y_mm * GeoModelKernelUnits::mm : casingSize);
  if (pitchX < casingSize || pitchY < casingSize)
    throw std::runtime_error("Module pitch is too small for the configured air and steel casing");

  const double tolX = std::max(cfg.tol_x_mm, std::abs(cfg.detector_offset_x_mm)) * GeoModelKernelUnits::mm;
  const double tolY = std::max(cfg.tol_y_mm, std::abs(cfg.detector_offset_y_mm)) * GeoModelKernelUnits::mm;
  const double worldXY = std::max((nx - 1) * pitchX + casingSize + 2.0 * tolX,
                                  (ny - 1) * pitchY + casingSize + 2.0 * tolY);
  const double worldZ = 120827 * GeoModelKernelUnits::mm;

  auto* worldShape = new GeoBox(0.5*worldXY, 0.5*worldXY, 0.5*worldZ);
  auto* worldLog   = new GeoLogVol("WorldLog", worldShape, air);
  auto* worldPhys  = new GeoPhysVol(worldLog);

  const double x0 = -0.5 * (nx - 1) * pitchX;
  const double y0 = -0.5 * (ny - 1) * pitchY;
  const double stackZ = aworldZ * GeoModelKernelUnits::mm;
  const double halfStackZ = 0.5 * stackZ;
  const double moduleWorldZ = cfg.detector_offset_z_mm * GeoModelKernelUnits::mm - 0.5 * worldZ;

  auto addCasingPart = [&](const std::string& name, double halfX, double halfY,
                           double x, double y, double z, GeoMaterial* material) {
    auto* shape = new GeoBox(halfX, halfY, halfStackZ);
    auto* log = new GeoLogVol((name + "_LOG").c_str(), shape, material);
    auto* phys = new GeoPhysVol(log);
    worldPhys->add(new GeoNameTag(name.c_str()));
    worldPhys->add(new GeoTransform(GeoTrf::Translate3D(x, y, z)));
    worldPhys->add(phys);
  };

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      const int mx = ix + 1;
      const int my = iy + 1;
      const double x = x0 + ix * pitchX + cfg.detector_offset_x_mm * GeoModelKernelUnits::mm;
      const double y = y0 + iy * pitchY + cfg.detector_offset_y_mm * GeoModelKernelUnits::mm;
      const double z = moduleWorldZ + (cfg.center_stack ? 0.0 : halfStackZ);
      const std::string moduleName = "MODULE_MX" + std::to_string(mx) + "Y" + std::to_string(my);

      auto* modShape = new GeoBox(0.5*plateSize, 0.5*plateSize, halfStackZ);
      auto* modLog   = new GeoLogVol((moduleName + "_LOG").c_str(), modShape, MM.air());
      auto* modPhys  = new GeoPhysVol(modLog);
      worldPhys->add(new GeoNameTag(moduleName.c_str()));
      worldPhys->add(new GeoTransform(GeoTrf::Translate3D(x, y, moduleWorldZ)));
      worldPhys->add(modPhys);
      auto moduleCfg = cfg;
      moduleCfg.module_nx = 1;
      moduleCfg.module_ny = 1;
      CalorimeterBuilder::buildStack(modPhys, MM, moduleCfg, mx, my);

      const double airCenter = 0.5 * plateSize + 0.5 * airThickness;
      const double steelCenter = 0.5 * plateSize + airThickness + 0.5 * steelThickness;
      const double airOuter = plateSize + 2.0 * airThickness;
      const double steelOuter = casingSize;
      struct CasingSide {
        const char* name;
        double dx;
        double dy;
        double airHalfX;
        double airHalfY;
        double steelHalfX;
        double steelHalfY;
      };
      const CasingSide sides[] = {
        {"Left",   -1.0,  0.0, 0.5 * airThickness, 0.5 * plateSize, 0.5 * steelThickness, 0.5 * airOuter},
        {"Right",   1.0,  0.0, 0.5 * airThickness, 0.5 * plateSize, 0.5 * steelThickness, 0.5 * airOuter},
        {"Bottom",  0.0, -1.0, 0.5 * airOuter, 0.5 * airThickness, 0.5 * steelOuter, 0.5 * steelThickness},
        {"Top",     0.0,  1.0, 0.5 * airOuter, 0.5 * airThickness, 0.5 * steelOuter, 0.5 * steelThickness}
      };
      for (const auto& side : sides) {
        addCasingPart(moduleName + "_Casing_Air_" + side.name,
                      side.airHalfX, side.airHalfY,
                      x + side.dx * airCenter, y + side.dy * airCenter, z, MM.air());
        addCasingPart(moduleName + "_Casing_Steel_" + side.name,
                      side.steelHalfX, side.steelHalfY,
                      x + side.dx * steelCenter, y + side.dy * steelCenter, z, MM.steel());
      }
    }
  }

  G4cout << "[DetectorConstruction] GeoModel world children = "
         << worldPhys->getNChildVols() << G4endl;

  return worldPhys;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Build GeoModel world (contains stack as children)
  m_geoWorld = buildGeoModelWorld();

  //  Create a Geant4 world (mother volume) 
  //  Convert GeoModel tree into an Assembly and imprint into lvWorld 
  bool ok = true;                 // builder doesn’t take ok in this API
  PVConstLink worldLink = m_geoWorld;
  
  // This builder is concrete (VolumeBuilder is abstract)
VolumeBuilder vb("World");   // key string can be anything; "World" is fine
  
  G4LogicalVolume* g4WorldLV = vb.Build(worldLink);
  
  if (!g4WorldLV) {
    G4Exception("DetectorConstruction::Construct",
                "Geo2G4BuildFailed",
                FatalException,
                "VolumeBuilder::Build returned null LV.");
  }
  
  auto* pvWorld = new G4PVPlacement(nullptr, {}, g4WorldLV, "WorldPV", nullptr, false, 0);
  
  G4cout << "[DetectorConstruction] After Build: world daughters = "
         << g4WorldLV->GetNoDaughters()
         << " LVStore size = " << G4LogicalVolumeStore::GetInstance()->size()
         << G4endl;
  
  

  //  Sensitive detector 
  auto* sdman = G4SDManager::GetSDMpointer();
  auto* caloSD = new CaloSD("CaloSD", m_store);
  sdman->AddNewDetector(caloSD);

  int nSensitive = 0;
  for (auto* lv : *G4LogicalVolumeStore::GetInstance()) {
    const std::string ln = lv->GetName();

    // Update these patterns to match your GeoLogVol names (print a few if needed)
    if (ln.find("Wide") != std::string::npos && ln.find("PVT") != std::string::npos) {
      lv->SetSensitiveDetector(caloSD); ++nSensitive;
    } else if (ln.find("Thin") != std::string::npos &&
              (ln.find("PS") != std::string::npos || ln.find("Poly") != std::string::npos)) {
      lv->SetSensitiveDetector(caloSD); ++nSensitive;
    } else if (ln.find("HPL") != std::string::npos &&
     ln.find("FiberCore") != std::string::npos) {
     lv->SetSensitiveDetector(caloSD); ++nSensitive;
    }
  }

  G4cout << "[DetectorConstruction] Sensitive LVs set: " << nSensitive << G4endl;
  auto* store = G4LogicalVolumeStore::GetInstance();
  
  const int visMode = m_vis_mode;

  auto* visHide = new G4VisAttributes(false);
  auto* visWorld = new G4VisAttributes(G4Colour(1, 1, 1, 0.05));
  visWorld->SetForceWireframe(true); visWorld->SetVisibility(true);
  auto* visModule = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.15));
  visModule->SetVisibility(true); visModule->SetForceWireframe(true);

  if (visMode == 0) {
    visModule->SetDaughtersInvisible(true);
    for (auto* lv : *store) {
      const auto& n = lv->GetName();
      if      (n == "WorldLog" || n == "World")          lv->SetVisAttributes(visWorld);
      else if (n.find("ModuleLog") != std::string::npos) lv->SetVisAttributes(visModule);
      else                                               lv->SetVisAttributes(visHide);
    }
    G4cout << "[DetectorConstruction] vis_mode=0: module envelopes only." << G4endl;

  }


  else if (visMode == 1) {
    auto* visECALlead  = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.7)); visECALlead->SetForceSolid(true);  visECALlead->SetDaughtersInvisible(true);
    auto* visECALscint = new G4VisAttributes(G4Colour(0.0, 0.8, 0.4, 0.7)); visECALscint->SetForceSolid(true); visECALscint->SetDaughtersInvisible(true);
    auto* visECALhpl   = new G4VisAttributes(G4Colour(0.1, 0.3, 1.0, 0.7)); visECALhpl->SetForceSolid(true);   visECALhpl->SetDaughtersInvisible(true);
    auto* visHCALiron  = new G4VisAttributes(G4Colour(0.9, 0.1, 0.1, 0.7)); visHCALiron->SetForceSolid(true);  visHCALiron->SetDaughtersInvisible(true);
    auto* visHCALscint = new G4VisAttributes(G4Colour(1.0, 0.6, 0.0, 0.7)); visHCALscint->SetForceSolid(true); visHCALscint->SetDaughtersInvisible(true);
    auto* visHCALhpl   = new G4VisAttributes(G4Colour(0.6, 0.0, 0.9, 0.7)); visHCALhpl->SetForceSolid(true);   visHCALhpl->SetDaughtersInvisible(true);
    auto* visCasingAir = new G4VisAttributes(G4Colour(0.2, 0.8, 1.0, 0.2)); visCasingAir->SetForceSolid(true);
    auto* visCasingSteel = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75, 1.0)); visCasingSteel->SetForceSolid(true);
    for (auto* lv : *store) {
      const auto& n = lv->GetName();
      if      (n == "WorldLog" || n == "World")          lv->SetVisAttributes(visWorld);
      else if (n.find("_Casing_Air_") != std::string::npos) lv->SetVisAttributes(visCasingAir);
      else if (n.find("_Casing_Steel_") != std::string::npos) lv->SetVisAttributes(visCasingSteel);
      else if (n.find("ModuleLog") != std::string::npos) lv->SetVisAttributes(visModule);
      else if (n.find("ECAL") != std::string::npos && n.find("_LOG") != std::string::npos) {
        if      (n.find("Lead") != std::string::npos) lv->SetVisAttributes(visECALlead);
        else if (n.find("HPL")  != std::string::npos) lv->SetVisAttributes(visECALhpl);
        else                                           lv->SetVisAttributes(visECALscint);
      }
      else if (n.find("HCAL") != std::string::npos && n.find("_LOG") != std::string::npos) {
        if      (n.find("Iron") != std::string::npos) lv->SetVisAttributes(visHCALiron);
        else if (n.find("HPL")  != std::string::npos) lv->SetVisAttributes(visHCALhpl);
        else                                           lv->SetVisAttributes(visHCALscint);
      }
      else {
        lv->SetVisAttributes(visHide);
        static int dbgCount = 0;
      }
    }
  } else if(visMode == 2) {
    auto* visWide   = MakeVis(0.0, 0.7, 0.7, 1.0);
    auto* visThin   = MakeVis(1.0, 0.55, 0.0, 1.0);
    auto* visFibre  = MakeVis(0.1, 0.2, 1.0, 1.0);
    auto* visLead   = MakeVis(0.5, 0.5, 0.5, 1.0);
    auto* visIron   = MakeVis(0.9, 0.0, 0.0, 1.0);
    auto* visAlCase = MakeVis(0.8, 0.8, 0.8, 0.15);
    for (auto* lv : *store) {
      const auto& n = lv->GetName();
      if      (n.find("WidePVT") != std::string::npos) lv->SetVisAttributes(visWide);
      else if (n.find("ThinPS")  != std::string::npos) lv->SetVisAttributes(visThin);
      else if (n.find("Fibre")   != std::string::npos) lv->SetVisAttributes(visFibre);
      else if (n.find("Lead") != std::string::npos || n.find("Pb") != std::string::npos) lv->SetVisAttributes(visLead);
      else if (n.find("Iron") != std::string::npos || n.find("Fe") != std::string::npos) lv->SetVisAttributes(visIron);
      else if (n.find("HPL")  != std::string::npos && n.find("Al") != std::string::npos) lv->SetVisAttributes(visAlCase);
    }
    G4cout << "[DetectorConstruction] vis_mode=2: full detail." << G4endl;
  }
  static G4Mutex gdmlMutex = G4MUTEX_INITIALIZER;
  static bool gdmlWritten = false;

  const char* doGdml = std::getenv("G4_EXPORT_GDML");

    // In MT jobs, ensure only one thread writes
    G4AutoLock lock(&gdmlMutex);
  if(write_gdml){
    if (!gdmlWritten && G4Threading::IsMasterThread()) {
      const char* out = std::getenv("G4_GDML_OUT");
      std::string outName = out ? out : "geometry.gdml";

      G4GDMLParser parser;
      parser.SetStripFlag(false); // keep names
      // 3rd argument: write schema location attributes
      parser.Write(outName, pvWorld, /*storeReferences=*/true);

      G4cout << "[DetectorConstruction] Wrote GDML to: " << outName << G4endl;
      gdmlWritten = true;
    }
  }


  return pvWorld;
}


// helper
static G4VisAttributes* MakeVis(double r, double g, double b, double a=1.0, bool solid=true) {
    auto* vis = new G4VisAttributes(G4Colour(r, g, b, a));
    vis->SetVisibility(true);
    vis->SetForceSolid(solid);
    return vis;
}
