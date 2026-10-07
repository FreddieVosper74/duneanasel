
#ifndef NDLARGEOMETRY_H
#define NDLARGEOMETRY_H

#include <vector>
#include <cmath>

namespace ana{

namespace NDLArGeo {

    // Structure defining the bounding box of an individual TPC volume
    struct Box {
        double x_min, x_max;
        double y_min, y_max;
        double z_min, z_max;

        bool Contains(double x, double y, double z) const {
            return (x >= x_min && x <= x_max &&
                    y >= y_min && y <= y_max &&
                    z >= z_min && z <= z_max);
        }
    };

    //define the outer boundaries of the detector volume (in cm)
    const double NDLArXLo = -346.9;
    const double NDLArXHi = 346.9;
    const double NDLArYLo = -215.5;
    const double NDLArYHi = 81.7;
    const double NDLArZLo = 418.2;
    const double NDLArZHi = 913.3;

    // define the full TPC Dimensions from YAML file
    const double TPC_DX = 46.788;
    const double TPC_DY = 297.6;
    const double TPC_DZ = 95.232;

    const double HALF_DX = TPC_DX / 2.0;
    const double HALF_DY = TPC_DY / 2.0;
    const double HALF_DZ = TPC_DZ / 2.0;

    // Center positions of all 70 TPCs [X, Y, Z]
    inline const std::vector<std::vector<double>> TPC_CENTERS = {
        {-276.2885, -66.8713, 465.7559}, {-323.7115, -66.8713, 465.7559},
        {-276.2885, -66.8713, 565.7559}, {-323.7115, -66.8713, 565.7559},
        {-276.2885, -66.8713, 665.7559}, {-323.7115, -66.8713, 665.7559},
        {-276.2885, -66.8713, 765.7559}, {-323.7115, -66.8713, 765.7559},
        {-276.2885, -66.8713, 865.7559}, {-323.7115, -66.8713, 865.7559},
        {-176.2885, -66.8713, 465.7559}, {-223.7115, -66.8713, 465.7559},
        {-176.2885, -66.8713, 565.7559}, {-223.7115, -66.8713, 565.7559},
        {-176.2885, -66.8713, 665.7559}, {-223.7115, -66.8713, 665.7559},
        {-176.2885, -66.8713, 765.7559}, {-223.7115, -66.8713, 765.7559},
        {-176.2885, -66.8713, 865.7559}, {-223.7115, -66.8713, 865.7559},
        { -76.2885, -66.8713, 465.7559}, {-123.7115, -66.8713, 465.7559},
        { -76.2885, -66.8713, 565.7559}, {-123.7115, -66.8713, 565.7559},
        { -76.2885, -66.8713, 665.7559}, {-123.7115, -66.8713, 665.7559},
        { -76.2885, -66.8713, 765.7559}, {-123.7115, -66.8713, 765.7559},
        { -76.2885, -66.8713, 865.7559}, {-123.7115, -66.8713, 865.7559},
        {  23.7115, -66.8713, 465.7559}, { -23.7115, -66.8713, 465.7559},
        {  23.7115, -66.8713, 565.7559}, { -23.7115, -66.8713, 565.7559},
        {  23.7115, -66.8713, 665.7559}, { -23.7115, -66.8713, 665.7559},
        {  23.7115, -66.8713, 765.7559}, { -23.7115, -66.8713, 765.7559},
        {  23.7115, -66.8713, 865.7559}, { -23.7115, -66.8713, 865.7559},
        { 123.7115, -66.8713, 465.7559}, {  76.2885, -66.8713, 465.7559},
        { 123.7115, -66.8713, 565.7559}, {  76.2885, -66.8713, 565.7559},
        { 123.7115, -66.8713, 665.7559}, {  76.2885, -66.8713, 665.7559},
        { 123.7115, -66.8713, 765.7559}, {  76.2885, -66.8713, 765.7559},
        { 123.7115, -66.8713, 865.7559}, {  76.2885, -66.8713, 865.7559},
        { 223.7115, -66.8713, 465.7559}, { 176.2885, -66.8713, 465.7559},
        { 223.7115, -66.8713, 565.7559}, { 176.2885, -66.8713, 565.7559},
        { 223.7115, -66.8713, 665.7559}, { 176.2885, -66.8713, 665.7559},
        { 223.7115, -66.8713, 765.7559}, { 176.2885, -66.8713, 765.7559},
        { 223.7115, -66.8713, 865.7559}, { 176.2885, -66.8713, 865.7559},
        { 323.7115, -66.8713, 465.7559}, { 276.2885, -66.8713, 465.7559},
        { 323.7115, -66.8713, 565.7559}, { 276.2885, -66.8713, 565.7559},
        { 323.7115, -66.8713, 665.7559}, { 276.2885, -66.8713, 665.7559},
        { 323.7115, -66.8713, 765.7559}, { 276.2885, -66.8713, 765.7559},
        { 323.7115, -66.8713, 865.7559}, { 276.2885, -66.8713, 865.7559}
    };

    // function to construct TPC boxes and return a vector containing the 70 boxes
    inline std::vector<Box> GetTPCBoxes(double margin_x = 0.0, double margin_y = 0.0, double margin_z = 0.0) {
        std::vector<Box> boxes;
        boxes.reserve(TPC_CENTERS.size());

        for (const auto& center : TPC_CENTERS) { //this for loop access the TPC_CENTERS vector
            Box b;
            b.x_min = center[0] - HALF_DX + margin_x;
            b.x_max = center[0] + HALF_DX - margin_x;
            b.y_min = center[1] - HALF_DY + margin_y;
            b.y_max = center[1] + HALF_DY - margin_y;
            b.z_min = center[2] - HALF_DZ + margin_z;
            b.z_max = center[2] + HALF_DZ - margin_z;
            boxes.push_back(b);
        }
        return boxes;
    }

    //construct the outer active volume as a box with the specified margins
    inline Box GetOuterBox(double margin_x = 0.0, double margin_y = 0.0, double margin_z = 0.0) {
        Box outer_box;
        outer_box.x_min = NDLArXLo + margin_x;
        outer_box.x_max = NDLArXHi - margin_x;
        outer_box.y_min = NDLArYLo + margin_y;
        outer_box.y_max = NDLArYHi - margin_y;
        outer_box.z_min = NDLArZLo + margin_z;
        outer_box.z_max = NDLArZHi - margin_z - 100;
        return outer_box;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Position Checking functions:

    // Checks if position (x, y, z) is inside the outer active volume box
    inline bool IsInOuterFiducialVolume(double x, double y, double z, 
                                        double margin_x = 0.0, double margin_y=0.0, double margin_z = 0.0){
    
    return GetOuterBox(margin_x, margin_y, margin_z).Contains(x, y, z);
    }

    // Checks if position (x, y, z) is inside any of the individual TPC boxes
    inline bool IsInBoxesFiducialVolume(double x, double y, double z, double margin_x = 0.0, double margin_y = 0.0, double margin_z = 0.0) {
        const std::vector<Box> boxes = GetTPCBoxes(margin_x, margin_y, margin_z);
        for (const auto& box : boxes) {
            if (box.Contains(x, y, z)) { // here you want to insert a "&& AV conatins"
                return true;
            }
        }
        return false;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Function which takes the following inputs:
    //vertex position x, y, z
    //TPC box FV margins x, y, z
    //Outer FV margins x, y, z

    inline bool IsInCombinedNDLArFV(double x, double y, double z, 
                                    const std::vector<double>& box_margins,
                                    const std::vector<double>& outer_margins){
        double o_mx = outer_margins.size() >0 ? outer_margins[0] : 0.0;
        double o_my = outer_margins.size() >1 ? outer_margins[1] : 0.0;
        double o_mz = outer_margins.size() >2 ? outer_margins[2] : 0.0;

        if (!IsInOuterFiducialVolume(x, y, z, o_mx, o_my, o_mz)) return false;

        double b_mx = box_margins.size() > 0? box_margins[0] : 0.0;
        double b_my = box_margins.size() > 1? box_margins[1] : 0.0;
        double b_mz = box_margins.size() > 2? box_margins[2] : 0.0;

        return IsInBoxesFiducialVolume(x, y, z, b_mx, b_my, b_mz);

    }

}
}
#endif // NDLARGEOMETRY_H
