///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This header defines the boolean-returning functions to be used in the Nu_muCC0Pi reco interaction selection for ND-LAr

//It contains functions for cut flow LVL2, 3, 4 . LVL1 (Fiducial volume) function can be found in NDLArGeometry.h

//This header also contains functions acting on truth interactions, which are then used in a Var definition
//The Var can be called to identify signal and background sources in the Nu_muCC0Pi selection, via truth-matching


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef NU_MU_CC0PI_FUNCTIONS_H
#define NU_MU_CC0PI_FUNCTIONS_H

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana{

namespace NDLArFunctions {

//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//reco selection functions from level 2 to level 4
// you can find level 1 (FV) cut in the geometry header

//reco selection level 2 criteria - intercation produces produces muon
bool ExactlyOnePrimaryMuon(const caf::SRInteractionProxy* sr)
{
    int nMuons = 0;

    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && p.pdg == 13)
            ++nMuons;
    }

    return nMuons == 1;
}


// reco election level 3 criteria - produced muon has minimum track 20cm
bool PrimMuMinTrackLength(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
      if (p.primary == 1 && p.pdg == 13 )
      {
        //define track length in here
          double trackstartX = p.start.x;
          double trackstartY = p.start.y;
          double trackstartZ = p.start.z;

          double trackendX = p.end.x;
          double trackendY = p.end.y;
          double trackendZ = p.end.z; 

          double tracklength = std::sqrt(std::pow((trackendX - trackstartX),2) + std::pow((trackendY - trackstartY),2) + std::pow((trackendZ - trackstartZ),2));
        
          if (tracklength >= 20.0){
            return true;
          }
          else{
            return false;
          }
      }
       
    }
    return false;

}

// reco selection level 4 criteria 

bool ProducesZeroPions(const caf::SRInteractionProxy* sr)
{
    int nPions = 0;

    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && p.pdg == 211)
            ++nPions;
    }

    return nPions == 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TRUTH functions, used for identifying true signal and background contaminants in the selection
//These functions are implamented in the Var definition below


//check true FV position
bool InFVTrue (const caf::SRTrueInteractionProxy* sr){

  const std::vector<double> kDefaultBoxMargins = {5.0, 5.0, 5.0}; //in cm
  const std::vector<double> kDefaultAVMargins = {25.0, 25.0, 25.0}; // in cm

    double x = sr->vtx.x;
    double y = sr->vtx.y;
    double z = sr->vtx.z;

  return NDLArGeo::IsInCombinedNDLArFV(x, y, z, kDefaultBoxMargins, kDefaultAVMargins);

}

//check true muon track length
bool PrimaryTrueMuonHasMinimumTrackLength(const caf::SRTrueInteractionProxy* truth)
{
    for (const auto& p : truth->prim)
    {
        if (p.pdg == 13)
        {
            const double dx = p.end_pos.x - p.start_pos.x;
            const double dy = p.end_pos.y - p.start_pos.y;
            const double dz = p.end_pos.z - p.start_pos.z;

            const double tracklength = std::sqrt(dx*dx + dy*dy + dz*dz);

            return tracklength >= 20.0;
        }
    }

    return false;
}


//Returns truth matched interaction for a reco interaction input 
const caf::SRTrueInteractionProxy* GetTruthMatch(const caf::SRInteractionProxy* sr)
{
    if (sr->truth.empty())
        return nullptr;

    // Find best-matched true interaction
    size_t tidx = 0;
    float maxOverlap = 0.0;

    for (size_t i = 0; i < sr->truthOverlap.size(); ++i) {
        if (sr->truthOverlap[i] > maxOverlap) {
            maxOverlap = sr->truthOverlap[i];
            tidx = i;
        }
    }

    // Get ancestor SRProxy
    const caf::SRProxy* proxy = sr->Ancestor<caf::SRProxy>();

    // Return matched true interaction
    return &proxy->mc.nu[sr->truth[tidx]];
}

//enumerates the interaction type
enum TruthInteractionType {
    kNuMuCC0PiFV = 0,
    kNuMuCCnPiFV = 1,
    kNuMu_NuMuBarCCincOOFV = 2,
    kNCinc       = 3,
    kNuMuBarCCincFV = 4,
    kNuECCinc    = 5,
    kNuEBarCCinc = 6,
    kNoMatch     = 7,
    kUncatagorisable = 8
};








//generate Var enumerates to correct interaction type 
const Var GetTruthInteractionType([](const caf::SRInteractionProxy* sr) -> double { 
    const caf::SRTrueInteractionProxy* truthInt = GetTruthMatch(sr);
    if (!truthInt) 
        return static_cast<double>(kNoMatch);

    // ------------------
    //basic truth info
    //-------------------

    const int nupdg = truthInt->pdg;
    const bool nu_iscc = truthInt->iscc;

    const bool inFV =
            InFVTrue(truthInt);
    
    const bool zeroPi = (truthInt->npip == 0 && truthInt->npim == 0 && truthInt->npi0 == 0);

    const bool hasLongMuon = PrimaryTrueMuonHasMinimumTrackLength(truthInt);

    const bool nPi = (truthInt->npip > 0 || truthInt->npim > 0 || truthInt->npi0 > 0);


    //------------------------------
    //NC check.
    //------------------------------

    if (!nu_iscc) {
        return static_cast<double>(kNCinc);
    }

    //------------------------------
    // nu_mu CC check
    //------------------------------

    else if (nupdg == 14) {


            //nu_mu CC interaction outside the FV
            if(!inFV && hasLongMuon){
                return static_cast<double>(kNuMu_NuMuBarCCincOOFV);
            }
            // nu_mu CC 0Pi with a truth muon track >= 20 cm
            if(zeroPi && hasLongMuon){
                return static_cast<double>(kNuMuCC0PiFV);
            }
            //nu_mu CC with at least one pion
            else if(nPi){

                return static_cast<double>( kNuMuCCnPiFV);
            }

            else {
                //nu_mu CC in FV by not nPi
                //this catches 0Pi as well as anything unusual
                return static_cast<double>(kUncatagorisable);
            }

    }

    //------------------------------
    //anti-nu_mu cc
    //------------------------------

    else if (nupdg == -14) {

            if (!inFV && hasLongMuon) {
                return static_cast<double>(kNuMu_NuMuBarCCincOOFV);
            }
        
        //anti-nu_mu CC inside FV
        return static_cast<double >(kNuMuBarCCincFV);


    }

    //------------------------------
    //nu_e CC check
    //------------------------------


    else if (nupdg == 12) {
            return static_cast<double>(kNuECCinc);
    }   


    //------------------------------
    //anti-nu_e CC check
    //------------------------------

    else if(nupdg == -12){
            return static_cast<double>(kNuEBarCCinc);
    }

    //-----------------------------
    //if we get here, we have a CC interaction that is not nu_mu or nu
    //-----------------------------
    else{
        return static_cast<double>(kUncatagorisable);
    }

});


    




} // namespace NDLArFunctions
} // namespace ana
#endif // NU_MU_CC0PI_FUNCTIONS_H



