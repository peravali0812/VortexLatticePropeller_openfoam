/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Namespace
    None

Class
    openPropVLM2D

Description
    This is the horizontal-axis wind propeller array actuator line model class.
    It will set up an array of various kinds of propellers (currently blades 
    only) within a flow field.  The blade rotation rate is set or calculated
    based on a simple torque control model (not implemented yet), the blades
    are rotated at each time step, the propeller is yawed (not implemented yet),
    the blade forces are calculated, the actuator line body force information
    is passed back to the flow solver, and propeller information is written to
    files.

SourceFiles
    openPropVLM2D.C

\*---------------------------------------------------------------------------*/

#ifndef openPropVLM2D_H
#define openPropVLM2D_H

#include "HashPtrTable.H"
#include "IOdictionary.H"
#include "IFstream.H"
#include "OFstream.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace propellerModels
{

/*---------------------------------------------------------------------------*\
                           Class openPropVLM2D declaration
\*---------------------------------------------------------------------------*/

class openPropVLM2D
{

private:
    // Private Data

        //- Constants
            //- Runtime pointer.
            const Time& runTime_;

            //- Mesh pointer.
            const fvMesh& mesh_;

            //- Velocity field pointer.
            const volVectorField& U_;

            //- Degrees to radians conversion factor.
            const scalar degRad;

            //- Revolutions per minute to radians per second conversion factor.
            const scalar rpmRadSec;



        //- Current time step size.
        scalar dt;
 
        //- Current simulation time.
        word time;

        //- Body force field applied to fluid by propeller.
        volVectorField bodyForce;



        //- Turbine Array Level Data (all variables start with lower case letter).
            //- Body force vector field created by propeller array.
            //- List of names of propellers in array.
            List<word> propellerName;

            //- Number of propellers in array.
            int numPropellers;

            //- List of names of propeller types in array.
            DynamicList<word> propellerType;
 
            //- List of locations of bases of propellers in array relative to origin (m).
            DynamicList<vector> baseLocation;

            //- List of number of actuator line points on blades of propellers
            //  in array (set to uniform for now--here for future upgrades). 
            DynamicList<int> numBladePoints; 

            //- List of description of actuator line point distribution types
            //  for each propeller.
            DynamicList<word> pointDistType;

            //- List of body force normalization parameters for each propeller.
            DynamicList<scalar> epsilon;
            DynamicList<scalar> smearRadius; //(m)

            //- List of factor to scale sphere of action radius for each propeller.
            DynamicList<scalar> sphereRadiusScalar;

            //- List of tip/root loss correction type for each propeller.
            DynamicList<word> tipRootLossCorrType;
	    
	    //- Rotor rotation direction as viewed from upwind.  Options are
	    //  "cw" for clockwise and "ccw" for counter-clockwise.
	    DynamicList<word> rotationDir;
	    
            //- Initial or fixed rotor speed (rpm).  A positive value means
	    //  clockwise rotation for a clockwise rotating propeller (see rotationDir
	    //  above) or counter-clockwise rotation for a counter-clockwise
	    //  rotating propeller.
            DynamicList<scalar> rotSpeed;

            //- Initial blade 1 azimuth angle (degrees) (looking from upwind to 
	    //  downwind, a positive azimuth angle makes a clockwise movement if
	    //  this is a clockwise rotating propeller (see rotationDir above) or
	    //  or a counterclockwise movement if this is a counter-clockwise
	    //  rotating propeller).
            DynamicList<scalar> azimuth;

            //- Initial or fixed nacelle yaw angle.  Direction that the propeller
            //  is pointed in cardinal directions (i.e. 0 = north, 90 = east, 
            //  180 = south, 270 = west) (degrees).  This is converted to radians
            //  in the more standard mathematical convention of 0 degrees on the 
            //  + x axis and positive degrees in the counter-clockwise direction.
            DynamicList<scalar> nacYaw;

	    //- Last time when output written.
	    DynamicList<scalar> lastOutputTime;

    	    //- Last time step when output written.
	    DynamicList<label> outputIndex;

            //- Number of distinct propellers in array.
            int numPropellersDistinct;

            //- List of distinct names of propeller types in array.
            DynamicList<word> propellerTypeDistinct;

            //- ID label given to each distinct type of propeller in the array.
            DynamicList<label> propellerTypeID;



        //- Turbine Level Data (all variables start with a capital letter).

            //*** THE FOLLOWING VARIABLES MATCH FAST INPUT FILE ***
            //- Number of blades;
            DynamicList<int> NumBl;

            //- Distance from rotor apex to blade tip (m).
            DynamicList<scalar> TipRad;

            //- Distance from rotor apex to blade root (m).
            DynamicList<scalar> HubRad;

            //- Advance ratio
            DynamicList<scalar> AdvanceRatio;

	    //- Ship speed
	    DynamicList<scalar> VS;

            //- Desired Ct
            DynamicList<scalar> DesiredCT;

            //- Distance from nacelle yaw axis to teeter pin or rotor apex (m).
            DynamicList<scalar> OverHang;

            //- Vertical distance from tower-top to rotor shaft centerline (m).
            DynamicList<scalar> Baseline2Shft;

            //- Shaft tilt-up angle (degrees).
            DynamicList<scalar> ShftTilt;

            //Number of section points(no of panels)*****************new one
            DynamicList<scalar> sectionPoints;

            // Start point  of the panel**************************new one
            //DynamicList<scalar> startPt;

            //End point of the panel*****************************new one
            //DynamicList<scalar> endPt;

            //- Coning angle of blades (degrees) (one for each blade).
            DynamicList<List<scalar> > Rake;

            //*** END OF FAST INPUT FILE VARIABLES ***

            //- Yaw rate of nacelle (NOT IMPLEMENTED).
            DynamicList<scalar> YawRate;

            //- Rotor speed controller type (NOT IMPLEMENTED, just runs at 
	    //  specified rotSpeed).
            DynamicList<word> SpeedControllerType;

	    //- Nacelle yaw controller type (NOT IMPLEMENTED, just remains at 
	    //  specified nacYaw).
	    DynamicList<word> YawControllerType;

            //- List of airfoils that compose propeller blade;
            DynamicList<List<word> > AirfoilType;

            //- Lists of blade data for each propeller type.
                //- Overall blade data array.
                DynamicList<List<List<scalar> > > BladeData;

                //- Station along blade in which information is given (m).
                DynamicList<DynamicList<scalar> > BladeStation;

                //- Chord at this station (m).
                DynamicList<DynamicList<scalar> > BladeChord;

		//- Circulation distribution (have to change this to vector it depends on our case code)
	        //DynamicList<scalarRectangularMatrix>  Circulation;
                RectangularMatrix<scalar> Circulation;

                //- Drag coefficient distribution
                DynamicList<DynamicList<scalar> > DragCoeff;

               //Pitch distribution__________________________newone
               DynamicList<DynamicList<scalar> > BladePitch;

              //Mean aerodynamic chord (remove this after u know abt sum of scalar field)_________new
               DynamicList<scalar> meanADchord;

              //Slope of cl vs alpha curve of the aerofoil------------------new
               DynamicList<scalar> a0;

              //camberSlope (dy/dx)________________________newone
               DynamicList<List<scalar> > camberSlope;

        //- Important Actuator Line Geometry Data.
            //- List of processors IDs corresponding to the processor domain that each
            //  actuator point lies within (a list of which processor is in control of
            //  of each actuator point).
            DynamicList<List<List<label> > > controlProcNo;

            //- List of cell ID that contains a certain actuator line point on the
            //  processor in control of that point.
            DynamicList<List<List<label> > > minDisCellID;

            //- List of locations of the intersection of the tower axis and the shaft 
            //  centerline relative to the origin (m).
            DynamicList<vector> boss;

            //- List of locations of the rotor apex relative to the origin (m).
            DynamicList<vector> rotorApex;

            //- List of list of labels or cells within sphere of action.
            DynamicList<DynamicList<label> > sphereCells;

            //- Actuator element width.
            DynamicList<DynamicList<scalar> > db;

            //- Control point locations with respect to origin.
            DynamicList<List<List<vector> > > bladePoints;

            //- Vortex point locations with respect to origin.
            DynamicList<List<List<vector> > > bladeVortexPoints;


            //- Total actuator line points in array.
            int totBladePoints;

            //- Blade radius away from rotor apex.  Must take into account coning.
            DynamicList<List<List<scalar> > > bladeRadius;

            //- Blade radius away from rotor apex.  Must take into account coning.
            DynamicList<List<List<scalar> > > bladeVortexRadius;


	    //- An indicator of shaft direction.  The convention is that when viewed
            //  from upwind, the rotor turns clockwise for positive rotation angles,
            //  regardless of if it is an upwind or downwind propeller.  uvShaft is
            //  found by subtracting the rotor apex location from the tower shaft
            //  intersection point.  This vector switches direciton depending on 
            //  if the propeller is upwind or downwind, so this uvShaftDir multiplier
            //  makes the vector consistent no matter what kind of propeller.
	    DynamicList<scalar> uvShaftDir;

            //- Unit vector pointing along the rotor shaft (axis of blade rotation).
            DynamicList<vector> uvShaft;

	    //- Unit vector pointing along the tower (axis of yaw).
	    DynamicList<vector> uvZ;
            DynamicList<vector> uvY;

            //- Blade force at each actuator point.
            DynamicList<List<List<vector> > > bladeForce;

            //- Three vectors for each blade of each propeller that define the local
            //  blade-aligned coordinate system.  Vector 0 is along the blade pointed
            //  from root to tip, vector 1 is in the tangential direction (direction
            //  of blade rotation) where positive points in the direction opposite 
            //  rotation if the rotor turns clockwise as viewed from upstream, and 
            //  points orthogonal to vector 0 and 1 and points toward downstream (but 
            //  vector 2 is not perfectly aligned with downstream due to rotor coning  
            //  and nacelle tilt).
            DynamicList<List<List<vector> > > bladeAlignedVectors;

            //- Wind vector at each actuator point in blade-aligned coordinate system.
            DynamicList<List<List<vector> > > windVectors;

            //- Induced velocity vector at each control point in blade-aligned coordinate system.
            DynamicList<List<List<vector> > > USTAR;

            //- Total inflow velocity vector at each control point in blade-aligned coordinate system.
            DynamicList<List<List<vector> > > VSTAR;

            //- Total inflow velocity vector at each control point in blade-aligned coordinate system.
            DynamicList<List<List<scalar> > > VSTARmag;

            //- Horseshoe influence matrix on control points
            DynamicList<List<List<List<vector> > > > UHIFc;


            //- Induced pitch angle at control points
            DynamicList<List<List<scalar> > > TANBIc;

            //- Induced pitch angle at vortices
            DynamicList<List<List<scalar> > > TANBIv;

            //- Wrench induced velocities
            DynamicList<List<List<vector> > > UW;

	    //- Change in yaw each time step.
	    DynamicList<scalar> deltaNacYaw;

	    //- Change in azimuth each time step.
	    DynamicList<scalar> deltaAzimuth;



        //- Information critical to propeller performance that can be written to file
        //  every time step.

            //- Angle of attack at each actuator point.
            DynamicList<List<List<scalar> > > alpha;

            //- Wind magnitude (not including radial wind) at each actuator point.
            DynamicList<List<List<scalar> > > Vmag;

            //- Lift/density at each actuator point. 
            DynamicList<List<List<vector> > > lift;

            //- Drag/density at each actuator point. 
            DynamicList<List<List<vector> > > drag;

            //- Axial force/density at each actuator point (not pointed in blade-local
            //  axial, but rather along shaft). 
            DynamicList<List<List<scalar> > > axialForce;

            //- Tangential force/density at each actuator point.
            DynamicList<List<List<scalar> > > tangentialForce;

            //- Thrust/density on propeller.
            DynamicList<scalar> thrust;

            //- Total torque/density on propeller.
            DynamicList<scalar> torque;

            //- Power/density on propeller.
            DynamicList<scalar> power;


    
    // Private Member Functions
        
        //- Rotate the blades.
        void rotateBlades();
        
        //- Yaw the nacelle.
        void yawNacelle();
        
        //- Calculate the blade rotation rate (fixed rate based on user-specified
	//  initial rotor speed or calculated with simple variable speed 
	//  torque-controller).
        void computeRotSpeed();

	//- Calculate the nacelle yaw position (fixed based on user-specified
	//  initial yaw direction or calculated with simple yaw controller).
	void computeNacYaw();

        //- Find out which processor zone each actuator line point lies within.
        //  (Which processor is controlling each actuator line point?)
        void findControlProcNo();
        
	//- Compute the unit vector at all blade points in blade-aligned
	//  coordinates.  Then, compute the wind vector in these blade-aligned
        //  coordinates at each actuator point.
	void computeWindVectors();

        //- Compute blade forces.
        void computeBladeForce();
        
        //- Compute body forces.
        void computeBodyForce();

        //- Rotates a point about a rotation axis and rotation point by the specified
        //  angle in radians.
        vector rotatePoint(vector point, vector rotationPoint, vector axis, scalar angle);

        //- Perform interpolation.
        scalar interpolate(scalar xNew, DynamicList<scalar>& xOld, DynamicList<scalar>& yOld);
        scalar interpolateField(scalar xNew, DynamicList<scalar>& xOld, scalarRectangularMatrix& yOld);
        label  interpolate(scalar xNew, DynamicList<scalar>& xOld, DynamicList<label>& yOld);

        //- Print variables for debugging.
        void printDebug();

	//- Estimate the initial indiced velocity
        void estimateUSTAR();

        //- Compute the magnitude of the total inflow velocity
        void computeVSTARmag();

	//- Compute tan betaI
        void findTanBetaI();

        //- Compute normalize Circulation
        void normalizeCirculation();

	//- Compute induced velocities
	void inducedVelocity();

	//- Compute horseshoe influence
	void horseShoe();

	//- Compute Wrench influence function
	vector wrench(scalar Z, scalar& TANBIv, scalar& Rc, scalar& Rv);

	//- Iterate induced velocities
        void iterateInducedVelocities();
//- Matrix Functions
scalarRectangularMatrix FillMatrix(Foam::scalarRectangularMatrix& M, label n);
//scalarRectangularMatrix FillMatrix(Foam::scalarSquareMatrix& M, label n);
scalarRectangularMatrix InverseMatrix(Foam::scalarRectangularMatrix& M, label n);
scalarRectangularMatrix GetMatrix(Foam::scalarRectangularMatrix& M, label n);
//scalarSquareMatrix GetMatrix(Foam::scalarRectangularMatrix& M, label n);
scalarRectangularMatrix MultiplyMatrixField(Foam::scalarRectangularMatrix& S, scalarRectangularMatrix& f, label n, label m);
//scalarRectangularMatrix MultiplyMatrixField(Foam::scalarSquareMatrix& S, scalarRectangularMatrix& f, label n, label m);
scalarRectangularMatrix sum_rows(Foam::scalarRectangularMatrix& S, label n, label m);
scalarRectangularMatrix MultiplyField(Foam::scalarRectangularMatrix& M, scalarField& f, label m);
public:
	
    //- Constructor
    openPropVLM2D
    (
        const volVectorField& U
    );
    
    
    //- Destructor
    virtual ~openPropVLM2D()
    {}
    
    
    // Public Member Functions
        //- Update state of propeller.
        void update();

	//- Return force.
	volVectorField& force();
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace propellerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

