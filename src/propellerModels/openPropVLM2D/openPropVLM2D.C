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

\*---------------------------------------------------------------------------*/

#include "openPropVLM2D.H"
namespace Foam
{
namespace propellerModels
{

// * * * * * * * * * * * * * *  Constructor  * * * * * * * * * * * * * * * * //

openPropVLM2D::openPropVLM2D
(
    const volVectorField& U
)	      
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity field
    U_(U),

    // Set the degrees to radians convesion factor.
    degRad((Foam::constant::mathematical::pi)/180.0),

    // Set the revolutions/min to radians/s conversion factor.
    rpmRadSec(2.0*(Foam::constant::mathematical::pi)/60.0),

    // Set the time step size.
    dt(runTime_.deltaT().value()),

    // Set the current simulation time.
    time(runTime_.timeName()),

    // Initialize the body force.
    bodyForce
    (
        IOobject
        (
            "bodyForce",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimForce/dimVolume/dimDensity,vector::zero)
    )

{   
    // Define dictionary that defines the propeller array.
    IOdictionary propellerArrayProperties
    (
        IOobject
        (
            "propellerArrayProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    


    // Read in the propeller array properties dictionary.  This is the uppermost level dictionary
    // that describes where the propellers are, what kind they are, their initial state, and 
    // information about how the lifting line method is applied to each propeller.
    propellerName = propellerArrayProperties.toc();

    numPropellers = propellerName.size();

    forAll(propellerName,i)
    {
        propellerType.append(word(propellerArrayProperties.subDict(propellerName[i]).lookup("propellerType")));
        baseLocation.append(vector(propellerArrayProperties.subDict(propellerName[i]).lookup("baseLocation")));
        numBladePoints.append(int(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("numBladePoints"))));
        pointDistType.append(word(propellerArrayProperties.subDict(propellerName[i]).lookup("pointDistType")));
        epsilon.append(scalar(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("epsilon"))));
        smearRadius.append(scalar(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("smearRadius"))));
        sphereRadiusScalar.append(scalar(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("sphereRadiusScalar"))));
        tipRootLossCorrType.append(word(propellerArrayProperties.subDict(propellerName[i]).lookup("tipRootLossCorrType")));
        rotationDir.append(word(propellerArrayProperties.subDict(propellerName[i]).lookup("rotationDir")));
        rotSpeed.append(scalar(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("RotSpeed"))));
        azimuth.append(scalar(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("Azimuth"))));
        nacYaw.append(scalar(readScalar(propellerArrayProperties.subDict(propellerName[i]).lookup("NacYaw")))); 
	lastOutputTime.append(runTime_.startTime().value());
	outputIndex.append(0);
    }




    // Catalog the various types of propellers.  For example if three propellers are GE 1.5 and 
    // two propellers are Siemens 2.3 machines, then assign the GE 1.5 an ID of 0 and the Siemens
    // 2.3 an ID of 1.
    numPropellersDistinct = 1;
    {
        propellerTypeDistinct.append(propellerType[0]);
        forAll(propellerType,i)
        {
            bool flag = false;
            for(int j = 0; j < numPropellersDistinct; j++)
            {
                if(propellerType[i] == propellerTypeDistinct[j])
                {
                   flag = true;
                }
            }
            if(flag == false)
            {
                numPropellersDistinct++;
                propellerTypeDistinct.append(propellerType[i]);
            }
        }
    }
    forAll(propellerType,i)
    {
        for(int j = 0; j < numPropellersDistinct; j++)
        {
            if(propellerType[i] == propellerTypeDistinct[j])
            {
                propellerTypeID.append(j);
            }
        }
    }




    // For each distinct propeller, read in properties of that propeller from separate
    // dictionaries.

    for(int i = 0; i < numPropellersDistinct; i++)
    {
        // Declare the propellerProperties dictionary for the ith propeller.
        IOdictionary propellerProperties
        (
            IOobject
            (
                propellerTypeDistinct[i],
                runTime_.constant(),"propellerProperties",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read in the data.
        Info << "Read in the data!" << endl;
        NumBl.append(scalar(readScalar(propellerProperties.lookup("NumBl"))));
        TipRad.append(scalar(readScalar(propellerProperties.lookup("TipRad"))));
        HubRad.append(scalar(readScalar(propellerProperties.lookup("HubRad"))));
        VS.append(scalar(readScalar(propellerProperties.lookup("Vs"))));
        DesiredCT.append(scalar(readScalar(propellerProperties.lookup("CTPDES"))));
        OverHang.append(scalar(readScalar(propellerProperties.lookup("OverHang"))));
        Baseline2Shft.append(scalar(readScalar(propellerProperties.lookup("Baseline2Shft"))));
        ShftTilt.append(scalar(readScalar(propellerProperties.lookup("ShftTilt"))));
        Rake.append(propellerProperties.lookup("Rake"));
        YawRate.append(scalar(readScalar(propellerProperties.lookup("YawRate"))));
        SpeedControllerType.append(word(propellerProperties.lookup("SpeedControllerType")));
        YawControllerType.append(word(propellerProperties.lookup("YawControllerType")));
        sectionPoints.append(scalar(readScalar(propellerProperties.lookup("sectionPoints"))));
        BladeData.append(propellerProperties.lookup("BladeData"));

        DynamicList<scalar> station;
        DynamicList<scalar> chord;
	DynamicList<scalar> Cd;
        DynamicList<scalar> Pitch;
        forAll(BladeData[i], j)
        {
            station.append(BladeData[i][j][0]);
            chord.append(BladeData[i][j][1]);
            Cd.append(BladeData[i][j][2]);
            Pitch.append(BladeData[i][j][3]);
        }

        BladeStation.append(station);
        BladeChord.append(chord);
        DragCoeff.append(Cd);
        BladePitch.append(Pitch);

        station.clear();
        chord.clear();
        Cd.clear();
        BladePitch.clear();

        meanADchord.append(scalar(readScalar(propellerProperties.lookup("meanADchord"))));
        a0.append(scalar(readScalar(propellerProperties.lookup("a0"))));
        camberSlope.append(propellerProperties.lookup("camberSlope"));
     }
//-----------------------------------------------------------------------------------------------
	    // Convert nacelle yaw from cardinal directions (like on a compass)
	    // to the normal convention of 0 degrees on the + x axis with
	    // positive degrees in the counter-clockwise direction.
	    forAll(nacYaw, i)
	    {
		if (nacYaw[i] > 180.0)
		{
		    nacYaw[i] = nacYaw[i] - 180.0;
		}
		else
		{
		    nacYaw[i] = nacYaw[i] + 180.0;
		}
		nacYaw[i] = 90.0 - nacYaw[i];
		if (nacYaw[i] < 0.0)
		{
		    nacYaw[i] = nacYaw[i] + 360.0;
		}
	    }
	    



	    // Convert quantities in degrees into radians (dynamic lists
	    // have to be done in loops).
	    azimuth  = degRad * azimuth;
	    rotSpeed = rpmRadSec * rotSpeed;
	    nacYaw   = degRad * nacYaw;
	    ShftTilt = degRad * ShftTilt;

	    forAll(Rake,i)
	    {
		Rake[i] = degRad * Rake[i];
	    }

	    // Calculate boss shaft intersection and rotor apex locations. (The
	    // i-index is at the propeller array level for each propeller and the j-
	    // index is for each type of propeller--if all propellers are the same, j-
	    // is always 0.)  The rotor apex is not yet rotated for initial yaw;
	    // that is done below.
	    for(int i = 0; i < numPropellers; i++)
	    {
		int j = propellerTypeID[i];
		boss.append(baseLocation[i]);
		boss[i].z() = boss[i].z() + Baseline2Shft[j];
		rotorApex.append(boss[i]);
		rotorApex[i].x() = rotorApex[i].x() + ((OverHang[j]) * Foam::cos(ShftTilt[j]));
		rotorApex[i].z() = rotorApex[i].z() +  (OverHang[j]) * Foam::sin(ShftTilt[j]);
	    }

	    // Define the cells that can possibly be influenced by the force
	    // exerted each propeller.  In otherwords, define a sphere of cell IDs
	    // around each propeller that will be saved into memory so that the
	    // entire domain need not be passed through when applying the force 
	    // field.  (The i-index is at the propeller array level for each 
	    // propeller, the j-index is for each type of propeller--if all propellers
	    // are the same, j is always 0, and the k-index is at the individual
	    // blade level.)
	    for(int i = 0; i < numPropellers; i++)
	    {
		DynamicList<label> sphereCellsI;
		scalar sphereRadius = 0.0;
		int j = propellerTypeID[i];
		forAll(Rake[j],k)
		{
		    scalar sphereRadiusI = sphereRadiusScalar[i] * Foam::sqrt(Foam::sqr(OverHang[j] + TipRad[j]*Foam::sin(Rake[j][k])) + Foam::sqr(TipRad[j]*Foam::cos(Rake[j][k])));
		    Info << "sphereRadiusI = " << sphereRadiusI << endl;
		    if(sphereRadiusI > sphereRadius)
		    {
			sphereRadius = sphereRadiusI;
		    }
		} 
		forAll(U_.mesh().cells(),cellI)
		{
		    if (mag(U_.mesh().C()[cellI] - boss[i]) <= sphereRadius)
		    {
			sphereCellsI.append(cellI);
		    }
		}
		sphereCells.append(sphereCellsI);
		sphereCellsI.clear();
	    }
	    
	    // Create the actuator line points (not yet rotated for initial nacelle***************
	    // yaw or initial rotor azimuth. i-index is at array level, j-index is
	    // for the type of propeller, k-index is for each blade, and m-index is
	    // for each actuator point.  Also create other important vectors, and
	    // initialize the blade force, blade aligned coordinate system, and
	    // wind vectors to zero.
	    totBladePoints = 0;
	    for(int i = 0; i < numPropellers; i++)
	    {
		int j = propellerTypeID[i];

		// Define which way the shaft points to distinguish between
		// upwind and downwind propellers.
		uvShaftDir.append(OverHang[j]/mag(OverHang[j]));

		// Define the vector along the shaft pointing in the
		// direction of the wind.
		uvShaft.append(rotorApex[i] - boss[i]);
		uvShaft[i] = (uvShaft[i]/mag(uvShaft[i])) * -uvShaftDir[i];
  
		// Define the vector aligned with the vertical coordinate pointing from
		// the ground to the boss.
		uvZ.append(boss[i] - baseLocation[i]);
		uvZ[i] = uvZ[i]/mag(uvZ[i]);

		uvY.append( uvShaft[i] ^ uvZ[i]);
		uvY.append( uvY[i] / mag(uvY[i]));
		// Calculate the length of each propeller section.
		db.append(DynamicList<scalar>(0));
		if(pointDistType[i] == "uniform")
		{
		    scalar liftingLineLength = (TipRad[j]-HubRad[j])/numBladePoints[i];
		    for(int m = 0; m < numBladePoints[i]; m++)
		    {
			db[i].append(liftingLineLength);
		    }
		}
		// Add other point distribution types here, such as cosine, tanh.

		// Now calculate the control points for each blade
		// of each propeller.  All blades points will be calculated
		// at zero azimuth (blade pointing up), and then rotated to its correct
		// position before doing a global rotation to the initial azimuth of
		// the rotor.  Also calculate the radius of each point.
		// Do the same for the vertices
		bladePoints.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));
		bladeVortexPoints.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i]+1,vector::zero)));
		bladeRadius.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));
		bladeVortexRadius.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i]+1,0.0)));

		for(int k = 0; k < NumBl[j]; k++)
		{
		    int tip = 0;
		    vector root = rotorApex[i];
		    vector apex = rotorApex[i];
		    scalar beta = Rake[j][k] - ShftTilt[j];
		    root.x() = root.x() + HubRad[j]*Foam::sin(beta);
		    root.z() = root.z() + HubRad[j]*Foam::cos(beta);
		    scalar dist = HubRad[j];
		    for(int m = 0; m < numBladePoints[i]; m++)
		    {
		       bladeVortexPoints[i][k][m].x() = apex.x() + dist*Foam::sin(beta);
		       bladeVortexPoints[i][k][m].y() = apex.y();
		       bladeVortexPoints[i][k][m].z() = apex.z() + dist*Foam::cos(beta);
		       bladeVortexRadius[i][k][m] = dist;               
		       dist = dist + 0.5*db[i][k];
		       bladePoints[i][k][m].x() = root.x() + dist*Foam::sin(beta);
		       bladePoints[i][k][m].y() = root.y();
		       bladePoints[i][k][m].z() = root.z() + dist*Foam::cos(beta);
		       bladeRadius[i][k][m] = dist;
		       totBladePoints++;
		       dist = dist + 0.5*db[i][k];
		       tip = m+1;
		    }
                    // include outermost vortex location 
                    bladeVortexPoints[i][k][tip].x() = apex.x() + dist*Foam::sin(beta);
		    bladeVortexPoints[i][k][tip].y() = apex.y();
		    bladeVortexPoints[i][k][tip].z() = apex.z() + dist*Foam::cos(beta);
		    bladeVortexRadius[i][k][tip] = dist; 
		    // Apply rotation to get blades, other than blade 1, in the right
		    // place.
		    if (k > 0)
		    {
			for(int m = 0; m < numBladePoints[i]; m++)
			{
			    bladePoints[i][k][m] = rotatePoint(bladePoints[i][j][m], rotorApex[i], uvShaft[i], (360.0/NumBl[j])*k*degRad);
			    bladeVortexPoints[i][k][m] = rotatePoint(bladeVortexPoints[i][j][m], rotorApex[i], uvShaft[i], (360.0/NumBl[j])*k*degRad);
                            tip = m+1;
			}
                        // include outermost vortex location 
                        bladeVortexPoints[i][k][tip] = rotatePoint(bladeVortexPoints[i][j][tip], rotorApex[i], uvShaft[i], (360.0/NumBl[j])*k*degRad);

		    }
		}

//-------------------------------------------- Influence Matrix[Aij]------------------------------------------------
//initilize xc and xv
scalarField xc(sectionPoints[j],0);
scalarField xv(sectionPoints[j],0);
// use the matrix algorithm initilize
scalarRectangularMatrix A(sectionPoints[j],sectionPoints[j],0);
for (int ni = 1; ni < sectionPoints[j]; ni++)
{
xc[ni] = 0.5 * (1 - (Foam::cos(ni*Foam::constant::mathematical::pi/sectionPoints[j])));
for (int mj = 1; mj < sectionPoints[j]; mj++)
{
xv[mj] = 0.5 * (1 - (Foam::cos((mj-0.5)*Foam::constant::mathematical::pi/sectionPoints[j])));
A[ni][mj] = (1/(2*3.141))/(xv[mj] - xc[ni]);
}
}
//----------------------------------Compute Bladepitch angle (theta)---------------------------------------------------
scalarField PitchAngle(numBladePoints[i],0.0);
//for (int m = 0; m < numBladePoints[i]; m++)
forAll(numBladePoints, m)
{
PitchAngle[m] = (Foam::atan(BladePitch[m][j]/Foam::constant::mathematical::pi)); 
}
//---------------------------------Compute angle beta (compute Vo and (omega*r)) --------------------------------------
//omega*r
scalarField omega_r(numBladePoints[i],0.0);
scalarField V_0(numBladePoints[i],0.0);
scalarField beta_initial(numBladePoints[i],0.0);
forAll(numBladePoints, m)
{
omega_r[m] = (rotSpeed[i]*rpmRadSec*(BladeStation[m][j])*(TipRad[j]-HubRad[j]));
V_0[m] = Foam::sqrt(Foam::sqr(VS[j]) + Foam::sqr(omega_r[m]));
beta_initial[m] = Foam::acos(omega_r[m]/V_0[m]);
}

//-------------------------Compute solidity(sigma) (include no of blades,mean AD chord)---------------------------------
scalarField sigma(numBladePoints[i],0.0);
forAll(numBladePoints,m)
{
sigma[m] = (NumBl[j]*(meanADchord[j]))/(Foam::constant::mathematical::pi*(BladeStation[j][m])*(TipRad[j]-HubRad[j]));
} 

//---------------------------------------------------Compute phi_0-----------------------------------------------------
scalarField phi_0(numBladePoints[i],0.0);
forAll(numBladePoints,m)
{
phi_0[m] = (PitchAngle[m]-beta_initial[m])/(1-(8*BladeStation[m][j]*Foam::sin(beta_initial[m])/(sigma[m]*a0[j])));
}

//-------------------------------------Compute beta_i, alpha_1 and V*_initial------------------------------------------
scalarField beta_i(numBladePoints[i],0.0);
scalarField alpha_i(numBladePoints[i],0.0);
scalarField Vstar_i(numBladePoints[i],0.0);
forAll(numBladePoints,m)
{
beta_i[m] = beta_initial[m] + phi_0[m];
alpha_i[m] = (PitchAngle[m]-beta_i[m])*(Foam::constant::mathematical::pi*180);
Vstar_i[m] = V_0[m]*Foam::cos(phi_0[m]);
}

//------------------------------Compute the induced velocities for each section----------------------------------------
scalarRectangularMatrix w_i(sectionPoints[j],numBladePoints[i],0.0);
forAll(numBladePoints,m)
{
forAll(sectionPoints,n)
{
w_i[n][m] = (camberSlope[j][n] - alpha_i[m])*Vstar_i[m];
}
}

//------Compute the circulation distribution along the chord and the blade, add  it inplace of assumed one-------------
//- Using Inverse Matrix
scalarRectangularMatrix Ainf = FillMatrix(A,sectionPoints[j]);
scalarRectangularMatrix InvAinf = InverseMatrix(Ainf,sectionPoints[j]);
scalarRectangularMatrix InvA = GetMatrix(InvAinf,sectionPoints[j]);
scalarRectangularMatrix g = MultiplyMatrixField(InvA,w_i,sectionPoints[j],numBladePoints[i]);
scalarRectangularMatrix gamma_i = sum_rows(g,1,numBladePoints[i]);
//-------------------Actual Circulation---------------------------
scalarField chord(numBladePoints[i],0.0);
forAll(numBladePoints,m)
{
chord[m] = (BladeChord[m][j]*2*TipRad[j]);
}

Circulation=MultiplyField(gamma_i,chord,numBladePoints[i]);
//Circulation.append(Circulation);
//---------------------------------------------------------------------------------------------------------------------------------------
		// Define the size of the bladeForce array and set to zero.
		bladeForce.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));
	  
		// Define the size of the bladeAlignedVectors array and set to zero.
		bladeAlignedVectors.append(List<List<vector> >(NumBl[j],List<vector>(3,vector::zero)));

		// Define the windVectors array and set it to zero.
		windVectors.append(List<List<vector> >(NumBl[j],List<vector>(numBladePoints[i],vector::zero)));


		// Define the induced velocity array and set it to zero.
		USTAR.append(List<List<vector> >(NumBl[j],List<vector>(numBladePoints[i],vector::zero)));

                // Define the total inflow velocity array and set it to zero.
                VSTAR.append(List<List<vector> >(NumBl[j],List<vector>(numBladePoints[i],vector::zero)));

                //- Define the magnitude array of the total inflow velocity
                VSTARmag.append(List<List<scalar> >(NumBl[j],List<scalar>(numBladePoints[i],0.0)));

		// Define the horseshoe influence function array and set it to zero.
                UHIFc.append(List<List<List<vector> > >(NumBl[j],List<List<vector> >(numBladePoints[i],List<vector>(numBladePoints[i], vector::zero))));

                // Define the horseshoe influence function array and set it to zero.
                //UHIFv.append(List<List<vector> >(NumBl[j],List<vector>(numBladePoints[i],vector::zero)));

		// Define the size of induced pitch angle lists at controlpoints and set to zero.
		TANBIc.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

		// Define the size of induced pitch angle lists at vortices and set to zero.
		TANBIv.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i]+1,0.0)));

                // Define the Wrench induced velocity array and set it to zero.
                UW.append(List<List<vector> >(NumBl[j],List<vector>(numBladePoints[i]+1,vector::zero)));
		// Define the size of the deltaNacYaw and deltaAzimuth lists and set to zero.
		deltaNacYaw.append(0.0);
		deltaAzimuth.append(0.0);

		// Define the size of the angle of attack lists and set to zero.
		alpha.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

		// Define the size of the wind speed magnitude lists and set to zero.
		Vmag.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

		// Define the size of the lift/density lists and set to zero.
                lift.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));

		// Define the size of the drag/density lists and set to zero.
		drag.append(List<List<vector> >(NumBl[j], List<vector>(numBladePoints[i],vector::zero)));

		// Define the size of the axial force/density lists and set to zero.
		axialForce.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

		// Define the size of the tangential force/density lists and set to zero.
		tangentialForce.append(List<List<scalar> >(NumBl[j], List<scalar>(numBladePoints[i],0.0)));

		// Define the size of the thrust/density lists and set to zero.
		thrust.append(0.0);

		// Define the size of the torque/density lists and set to zero.
		torque.append(0.0);

		// Define the size of the power/density lists and set to zero.
		power.append(0.0);

		// Define the size of the processor-in-control list and set to zero.
		controlProcNo.append(List<List<label> >(NumBl[j], List<label>(numBladePoints[i],0)));

		// Define the size of the cell-containing-actuator-point ID list and set to zero.
		minDisCellID.append(List<List<label> >(NumBl[j], List<label>(numBladePoints[i],0)));
	    }

	    // Yaw the nacelle to initial position.
	    deltaNacYaw = nacYaw;
	    yawNacelle();

	    // Rotate the rotor to initial azimuth angle.
	    deltaAzimuth =  azimuth;
	    rotateBlades();

	    // Find out which processors control each actuator line point.
	    findControlProcNo();

	    // Compute the wind vectors at this initial time step.
	    computeWindVectors();

            iterateInducedVelocities();

            computeVSTARmag();
            // Circulation computation done! - Now, compute forces on the blade

	    // Compute the blade forces due to this wind at the initial time step.
	    computeBladeForce();

	    // Compute the resultant body force at this initial time step.
            computeBodyForce();

	}
// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

	void openPropVLM2D::rotateBlades()
	{  
	    // Perform rotation propeller by propeller.
	    forAll(uvShaft, i)
	    {	
		// Check the rotation direction first and set the local delta azimuth
		// variable accordingly.
		scalar deltaAzimuthI = 0.0;
		if (rotationDir[i] == "cw")
		{
		    deltaAzimuthI =  deltaAzimuth[i];
		}
		if (rotationDir[i] == "ccw")
		{
		    deltaAzimuthI = -deltaAzimuth[i];
		}

		// Rotate propeller blades, blade by blade, point by point.
		forAll(bladePoints[i], j)
		{
		    forAll(bladePoints[i][j], k)
		    {
			bladePoints[i][j][k] = rotatePoint(bladePoints[i][j][k], rotorApex[i], uvShaft[i], deltaAzimuthI);
		    }
		}   

		// Calculate the new azimuth angle and make sure it isn't
		// bigger than 2*pi.
		azimuth[i] = azimuth[i] + deltaAzimuth[i];
		if (azimuth[i] >= 2.0 * Foam::constant::mathematical::pi)
		{
		    azimuth[i] -= 2.0 * Foam::constant::mathematical::pi;
		}
	    }
	}
		

	void openPropVLM2D::yawNacelle()
	{
	    // Perform rotation propeller by propeller.
	    forAll(uvZ, i)
	    {
		// Rotate the rotor apex first.
		rotorApex[i] = rotatePoint(rotorApex[i], boss[i], uvZ[i], deltaNacYaw[i]);

		// Recompute the shaft unit vector since the shaft has rotated.
		uvShaft[i] = rotorApex[i] - boss[i];
                uvShaft[i] = (uvShaft[i]/mag(uvShaft[i])) * -uvShaftDir[i];
		// Rotate propeller blades, blade by blade, point by point.
		forAll(bladePoints[i], j)
		{
		    forAll(bladePoints[i][j], k)
		    {
			bladePoints[i][j][k] = rotatePoint(bladePoints[i][j][k], boss[i], uvZ[i], deltaNacYaw[i]);
		    }
		}   

		// Compute the new yaw angle and make sure it isn't
		// bigger than 2*pi.
		nacYaw[i] = nacYaw[i] + deltaNacYaw[i];
		if (nacYaw[i] >= 2.0 * Foam::constant::mathematical::pi)
		{
		    nacYaw[i] -= 2.0 * Foam::constant::mathematical::pi;
		}
	    }
	}


	void openPropVLM2D::findControlProcNo()
	{
	    // Create two lists of lists.  The first is a list for each processor
	    // of lists for each actuator line point that contains minimum distance
	    // between the actuator line point and a CFD cell center.  The second
	    // is a list for each processor of lists for each actuator line point
	    // that contains the cell ID of the minimum distance CFD cell.
	    List<List<scalar> > minDisLocal(Pstream::nProcs(), List<scalar>(totBladePoints,1.0E30));
	    List<List<label> > minDisCellIDLocal(Pstream::nProcs(), List<label>(totBladePoints,-1));

	    int iter = 0;
	    forAll(bladePoints, i)
	    {
		forAll(bladePoints[i], j)
        	{
		    forAll(bladePoints[i][j], k)
		    {
			if (sphereCells[i].size() > 0)
			{
			// Find the cell that the actuator line point lies within and the distance
			// from the actuator line point to that cell center.
			label minDisCellID = sphereCells[i][0];
			scalar minDis = mag(mesh_.C()[minDisCellID] - bladePoints[i][j][k]);

			forAll(sphereCells[i], m)
			{
			    scalar dis = mag(mesh_.C()[sphereCells[i][m]] - bladePoints[i][j][k]);
			    if(dis <= minDis)
			    {
				minDisCellID = sphereCells[i][m];
			    }
			    minDis = mag(mesh_.C()[minDisCellID] - bladePoints[i][j][k]);
			}
			minDisLocal[Pstream::myProcNo()][iter] = minDis;
			minDisCellIDLocal[Pstream::myProcNo()][iter] = minDisCellID;
			}
			iter++;
		    }
		}
	    }

	    // Perform a parallel gather of these local lists to the master processor and
	    // and then parallel scatter the lists back out to all the processors.
	    Pstream::gatherList(minDisLocal);
	    Pstream::scatterList(minDisLocal);
	    Pstream::gatherList(minDisCellIDLocal);
	    Pstream::scatterList(minDisCellIDLocal);

	    //for(int i = 0; i<totBladePoints; i++)
	    //{
	    //    for(int j = 0; j<Pstream::nProcs(); j++)
	    //    {
	    //        Info << minDisLocal[j][i] << " ";
	    //    }
	    //    Info << endl;
	    //}

	    // For each actuator line point, find the minDisLocal list entry that is smallest.
	    // The processor to which this entry belongs is the processor that controls this
	    // actuator line point.  Then extract the corresponding entry from the minDisCellIDLocal
	    // list and put into the global list.
	    iter = 0;
	    forAll(controlProcNo, i)
	    {
		forAll(controlProcNo[i], j)
		{
		    forAll(controlProcNo[i][j], k)
		    {
			scalar minDis = 1.0E30;
			for(int m = 0; m < Pstream::nProcs(); m++)
			{
			    if(minDisLocal[m][iter] <= minDis)
			    {
				minDis = minDisLocal[m][iter];
				controlProcNo[i][j][k] = m;
				minDisCellID[i][j][k] = minDisCellIDLocal[m][iter];
			    }
			}
			iter++;
		    }
		}
	    }
	}
		

	void openPropVLM2D::computeRotSpeed()
	{
	    // Proceed propeller by propeller.
	    forAll(rotSpeed, i)
	    {
		int j = propellerTypeID[i];
		if (SpeedControllerType[j] == "none")
		{
		    // Do nothing.
		    deltaAzimuth[i] = rotSpeed[i] * dt;
		}
		else if (SpeedControllerType[j] == "simple")
		{
		    // Placeholder for when this is implemented.
	}
    }
}


void openPropVLM2D::computeNacYaw()
{
    // Proceed propeller by propeller.
    forAll(deltaNacYaw, i)
    {
        int j = propellerTypeID[i];
        if (YawControllerType[j] == "none")
        {
            // Do nothing.
	    deltaNacYaw[i] = 0.0;
        }
        else if (YawControllerType[j] == "simple")
        {
            // Placeholder for when this is implemented.
        }
    }
}	


void openPropVLM2D::computeWindVectors()
{
    // Create a list for each processor of lists for each actuator line point 
    // that contains wind velocity in x, y, z coordinates.  This list will
    // be parallel gathered/scattered to get the correct wind at all locations.
    List<List<vector> > windVectorsLocal(Pstream::nProcs(), List<vector>(totBladePoints,vector::zero));

    int iter = 0;
    forAll(bladePoints, i)
    {
        forAll(bladePoints[i], j)
        {
            forAll(bladePoints[i][j], k)
            {
                if(Pstream::myProcNo() == controlProcNo[i][j][k])
                {
                    windVectorsLocal[Pstream::myProcNo()][iter] = U_[minDisCellID[i][j][k]];
                }
                iter++;
            }
        }
    }

    // Perform a parallel gather of this local list to the master processor and
    // and then parallel scatter the list back out to all the processors.
    Pstream::gatherList(windVectorsLocal);
    Pstream::scatterList(windVectorsLocal);

    // Proceed propeller by propeller.
    iter = 0;
    forAll(windVectors, i)
    {
        int n = propellerTypeID[i];
        // Proceed blade by blade.
        forAll(windVectors[i], j)
        {
            // If clockwise rotating, this vector points along the blade toward the tip.
	    // If counter-clockwise rotating, this vector points along the blade toward the root.
	    if (rotationDir[i] == "cw")
	    {
                bladeAlignedVectors[i][j][2] =   bladePoints[i][j][0] - rotorApex[i];
                bladeAlignedVectors[i][j][2] =   bladeAlignedVectors[i][j][2]/mag(bladeAlignedVectors[i][j][2]);
	    }
	    else if (rotationDir[i] == "ccw")
	    {
                bladeAlignedVectors[i][j][2] = -(bladePoints[i][j][0] - rotorApex[i]);
                bladeAlignedVectors[i][j][2] =   bladeAlignedVectors[i][j][2]/mag(bladeAlignedVectors[i][j][2]);
	    }

            bladeAlignedVectors[i][j][1] = bladeAlignedVectors[i][j][2]^uvShaft[i];
            bladeAlignedVectors[i][j][1] = bladeAlignedVectors[i][j][1]/mag(bladeAlignedVectors[i][j][1]);

            bladeAlignedVectors[i][j][0] = bladeAlignedVectors[i][j][1]^bladeAlignedVectors[i][j][2];
            bladeAlignedVectors[i][j][0] = bladeAlignedVectors[i][j][0]/mag(bladeAlignedVectors[i][j][0]);

            // Proceed point by point.
            forAll(windVectors[i][j], k)
            {
                // Zero the wind vector.
                windVectors[i][j][k] = vector::zero;

                // Now put the velocity in that cell into blade-oriented coordinates.
                windVectors[i][j][k].x() = (bladeAlignedVectors[i][j][0] & windVectorsLocal[controlProcNo[i][j][k]][iter]) / VS[i];
                windVectors[i][j][k].y() = (bladeAlignedVectors[i][j][1] & windVectorsLocal[controlProcNo[i][j][k]][iter]) / VS[i] + (rotSpeed[i] * bladeRadius[i][j][k] * cos(Rake[n][j]));
                windVectors[i][j][k].z() = (bladeAlignedVectors[i][j][2] & windVectorsLocal[controlProcNo[i][j][k]][iter]) / VS[i];
                iter++;
            }
            
        }
    }
}


void openPropVLM2D::computeBladeForce()
{
    // Proceed propeller by propeller.
    forAll(windVectors, i)
    {
        int m = propellerTypeID[i];

        // Set the total thrust of the propeller to zero.  Thrust will be summed on a blade-element-
        // wise basis.
        thrust[i] = 0.0;

        // Set the total thrust of the propeller to zero.  Thrust will be summed on a blade-element-
        // wise basis.
        torque[i] = 0.0;

        // Proceed blade by blade.
        forAll(windVectors[i], j)
        {
            scalar bladeThrust = 0.0;
            // Proceed point by point.
            forAll(windVectors[i][j], k)
            {
                // Interpolate the local chord.
                scalar localChord = interpolate(bladeRadius[i][j][k], BladeStation[m], BladeChord[m]);

		// Interpolate the local Circulation
		scalar Gamma = interpolateField(bladeRadius[i][j][k], BladeStation[m], Circulation);
                vector bladeUp = vector::zero;
                bladeUp.z() = 1;
		// Interpolate the local given drag coefficient
		scalar localDrag = interpolate(bladeRadius[i][j][k], BladeStation[m], DragCoeff[m]);

                // Get the angle of the wind with respect to rotor plane tangent direction.
                scalar windAng = Foam::atan2(VSTAR[i][j][k].x(),VSTAR[i][j][k].y())/degRad; 

                // Apply tip/root-loss correction factor.
                // Tip/root-loss correction factor of Glauert.
                scalar F = 1.0;

                if(tipRootLossCorrType[i] == "none")
                {
                    F = 1.0;
                }

                else if(tipRootLossCorrType[i] == "Glauert")
                {
                    scalar g = 1.0;

                    scalar ftip  = (TipRad[m] - bladeRadius[i][j][k])/(bladeRadius[i][j][k] * sin(windAng*degRad));
                    scalar Ftip  = (2.0/(Foam::constant::mathematical::pi)) * acos(exp(-g * (NumBl[m] / 2.0) * ftip));

                    scalar froot = (bladeRadius[i][j][k] - HubRad[i])/(bladeRadius[i][j][k] * sin(windAng*degRad));
                    scalar Froot = (2.0/(Foam::constant::mathematical::pi)) * acos(exp(-g * (NumBl[m] / 2.0) * froot));

                    F = Ftip * Froot;
                }

		// Using the actual circulation, chord and profile drag coefficient to calculate
		// lift and drag per density in local blade coordinates
		// The lift vector is perpendicular to both, the inflow V*-vector and the normal 
		// vector of the local blade-section-plane
                
		lift[i][j][k] = F * VSTAR[i][j][k] ^( Gamma * bladeUp);
		drag[i][j][k] = -VSTAR[i][j][k] * (0.5 * F * Foam::sqr(VSTARmag[i][j][k]) * localDrag * localChord);

                // Make the scalar lift and drag quantities vectors in the Cartesian coordinate system.
                vector trueLift;
                trueLift.x() = lift[i][j][k] & bladeAlignedVectors[i][j][0];
                trueLift.y() = lift[i][j][k] & bladeAlignedVectors[i][j][1];
                trueLift.z() = lift[i][j][k] & bladeAlignedVectors[i][j][2];

                vector trueDrag;
                trueDrag.x() = drag[i][j][k] & bladeAlignedVectors[i][j][0];
                trueDrag.y() = drag[i][j][k] & bladeAlignedVectors[i][j][1];
                trueDrag.z() = drag[i][j][k] & bladeAlignedVectors[i][j][2];

                // Add up lift and drag to get the resultant force/density applied to this blade element.
                bladeForce[i][j][k] = -(trueDrag + trueLift);

                // Find the component of the blade element force/density in the axial (along the shaft)
                // direction.
                axialForce[i][j][k] = -bladeForce[i][j][k] & uvShaft[i];

                // Find the component of the blade element force/density in the tangential (torque-creating)
                // direction.
                tangentialForce[i][j][k] = bladeForce[i][j][k] & bladeAlignedVectors[i][j][1];

                // Add this blade element's contribution to thrust to the total propeller thrust.
                thrust[i] += axialForce[i][j][k];
                bladeThrust += axialForce[i][j][k];
                // Add this blade element's contribution to torque to the total propeller torque.
                torque[i] += tangentialForce[i][j][k] * bladeRadius[i][j][k] * cos(Rake[m][j]);
            }

        }
        // Compute power based on torque and rotation speed.
        power[i] = torque[i] * rotSpeed[i];
    }
}


void openPropVLM2D::computeBodyForce()
{  
    // Proceed propeller by propeller.
    forAll(bladeForce, i)
    {
        // Proceed over all sphere cells of that propeller if there are sphere cells.
        if (sphereCells[i].size() > 0)
        {
            forAll(sphereCells[i], m)
            {
                // For each sphere cell, set the body force to zero initially (this probably should be changed in
                // the future since propellers could forseeably be put so close together that this would zero out
                // the effect of the neighboring propeller.
                bodyForce[sphereCells[i][m]] = vector::zero;

                // For each blade.
                forAll(bladeForce[i], j)
                {   
                    // For each blade point.
                    forAll(bladeForce[i][j], k)
                    {
                        scalar dis = mag(mesh_.C()[sphereCells[i][m]] - bladePoints[i][j][k]);
                        
                        if (dis <= smearRadius[i])
                        {

                            bodyForce[sphereCells[i][m]] += bladeForce[i][j][k] * (Foam::exp(-Foam::sqr(dis/epsilon[i]))/(Foam::pow(epsilon[i],3)*Foam::pow(Foam::constant::mathematical::pi,1.5)));

                        }
                    }
                }  
            }
        }
    }
}


vector openPropVLM2D::rotatePoint(vector point, vector rotationPoint, vector axis, scalar angle)
{
    // Declare and define the rotation matrix.
    tensor RM;
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle); 
    RM.xy() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle); 
    RM.xz() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle); 
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract rotation point
    // off the point to be rotated.
    point = point - rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation point.
    point = point + rotationPoint;

    return point;
}


scalar openPropVLM2D::interpolate(scalar xNew, DynamicList<scalar>& xOld, DynamicList<scalar>& yOld)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
            return yOld[indexM] - ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
            return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
            return yOld[indexP] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
            return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}
//--------------------------------------InterpolateField-----------------------------------------------------
scalar openPropVLM2D::interpolateField(scalar xNew, DynamicList<scalar>& xOld, scalarRectangularMatrix& yOld)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
            return yOld[1][indexM] - ((yOld[1][indexP] - yOld[1][indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
            return yOld[1][indexM] + ((yOld[1][indexP] - yOld[1][indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
    }
else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
            return yOld[1][indexP] + ((yOld[1][indexP] - yOld[1][indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
            return yOld[1][indexM] + ((yOld[1][indexP] - yOld[1][indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
    }
    else if (xNew == xOld[index])
    {
        return yOld[1][index];
    }
else
    {
        return 0.0;
    }
}



//-----------------------------------------------------------------------------------------------------------
label openPropVLM2D::interpolate(scalar xNew, DynamicList<scalar>& xOld, DynamicList<label>& yOld)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
            return yOld[indexM] - ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
            return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
            return yOld[indexP] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
            return yOld[indexM] + ((yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM]))*(xNew - xOld[indexM]);
        }
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}

//- Estimate the initial induced velocity
void openPropVLM2D::estimateUSTAR()
{	
    // for all propellers
    forAll(USTAR, i)
    {
        // for all blades
        forAll(USTAR[i], j)
        {       
            // for all panel on the lifting line
            forAll(USTAR[i][j], k)
            {
                USTAR[i][j][k].x() = 0.5 * (Foam::sqrt(1 + DesiredCT[i])-1) ;
            }
        }
    }
}
//- Compute the magnitude of the total inflow velocity
void openPropVLM2D::computeVSTARmag()
{ 
    forAll(VSTARmag, i)
    {
        forAll(VSTARmag[i], j)
        {
            forAll(VSTARmag[i][j], k)
            {
                VSTARmag[i][j][k] = Foam::sqrt(sqr(windVectors[i][j][k].x()+USTAR[i][j][k].x()) + sqr(windVectors[i][j][k].y() + USTAR[i][j][k].y())) ;
            }
        }
    }
}

//- Compute tan betaI
void openPropVLM2D::findTanBetaI()
{
    forAll(TANBIc, i)
    {
        forAll(TANBIc[i], j)
        {
            forAll(TANBIc[i][j], k)
            {
                VSTAR[i][j][k].x() = windVectors[i][j][k].x() + USTAR[i][j][k].x();
                VSTAR[i][j][k].y() = windVectors[i][j][k].y() + USTAR[i][j][k].y();

                TANBIc[i][j][k]    = VSTAR[i][j][k].x() / VSTAR[i][j][k].y() ;
            }
            int iter = 0;
            forAll(TANBIv[i][j], k)
	    // Assumption: 	The TANBI is the same at the outermost vortex points
	    //			as on the control points!
	    // LATER!		Implement an extrapolation from the control points
	    // 			to the vortex points!
            {
                DynamicList<scalar> tmp1;
                tmp1.append(bladeRadius[i][j]);
                DynamicList<scalar> tmp2;
                tmp2.append(TANBIc[i][j]);
                TANBIv[i][j][k]    = interpolate(bladeVortexRadius[i][j][k],tmp1, tmp2);
                tmp1.clear();
                tmp2.clear();
                iter++;
            }
            // Now, for the moment, set the outermost TANBIv to the outermost TANBIc
            TANBIv[i][j][0] = TANBIc[i][j][0];
            TANBIv[i][j][iter] = TANBIc[i][j][iter-1];
        }
    }
}
//- Compute induced velocities
void openPropVLM2D::inducedVelocity()
{   
    vector tmp = vector::zero;
    forAll(USTAR, i)
    {
        forAll(USTAR[i], j)
        {
            forAll(USTAR[i][j], k)
            {
                tmp = vector::zero;
                forAll(UHIFc[i][j][k], m)
                {
                    // Interpolate given circulation from input station to 
                    // bladePointRadius
                    scalar G = interpolateField(bladeRadius[i][j][m], BladeStation[i], Circulation) / (2 * Foam::constant::mathematical::pi * TipRad[i] * VS[i]);
                    tmp = tmp + UHIFc[i][j][k][m] * G ;

                }

                USTAR[i][j][k] =  tmp;
 
            }
        }
    }
}

//- Iterate induced velocities
void openPropVLM2D::iterateInducedVelocities()
{   
    int x = 0;
    estimateUSTAR();
    findTanBetaI();
    while ( x < 5)
    {
        findTanBetaI();
        // compute horseshoe influence 
        horseShoe();
        // compute VSTAR derivative NOT IMPLEMENTED

        // start circulation distribution optimisation
        // set up simultaneous equations for G and LM
        // and solve linear system of equation -NOT IMPLEMENTED

        // compute induced velocities
        inducedVelocity();

        // compute VSTAR magnitude
        computeVSTARmag();

        findTanBetaI();
        x++;
    }    
}

//- Compute horseshoe influence
void openPropVLM2D::horseShoe()
{
    forAll(UHIFc, i)
    {
        forAll(UHIFc[i], j)
        {
            forAll(UHIFc[i][j], k)
            {
                // for all vortices find induced velocity at RC by a 
                // unit vortex shed by RV(m) 
                // wrench returns 2*pi*R*u_bar
                forAll(UW[i][j], l)
                {
                    scalar normedRc = bladeRadius[i][j][k]/TipRad[i];
                    scalar normedRv = bladeVortexRadius[i][j][l]/TipRad[i];

                    UW[i][j][l] = wrench(NumBl[i],TANBIv[i][j][l], normedRc, normedRv);

                    // LATER! Include hub-image effects!!!
                    // LATER! Include swirl cancellation factor
                    // UW[i][j][k].y() = UW[i][j][k].y() * SCF
                }
                
                forAll(UHIFc[i][j][k], m)
                {
                    UHIFc[i][j][k][m] = UW[i][j][m+1] - UW[i][j][m];
                }
            }
        }
     }
 
}

//- Compute Wrench influence function
vector openPropVLM2D::wrench(scalar Z, scalar& TANBIv, scalar& Rc, scalar& Rv)
{
    scalar y  = Rc/(Rv*TANBIv);
    scalar y0 = 1/TANBIv;
    scalar U  = Foam::pow(((y0 * (Foam::sqrt(1+Foam::sqr(y))-1))*(Foam::exp(Foam::sqrt(1+Foam::sqr(y))-Foam::sqrt(1+ Foam::sqr(y0)))) / (y*(Foam::sqrt(1+Foam::sqr(y0))-1))),Z);
    vector uw = vector::zero; 
    if (Rc == Rv)
    {
        uw = vector::zero;
    }
    else if (Rc < Rv)
    {   
        scalar F1 =(-1 / (2*Z*y0)) * (Foam::pow(((1+Foam::sqr(y0))/(1+Foam::sqr(y))),0.25)) * ((U/(1-U)) + (1/(24*Z)*((9*Foam::sqr(y0)+2)/(Foam::pow((1+Foam::sqr(y0)),1.5)) + (3*Foam::sqr(y)-2)/(Foam::pow((1+Foam::sqr(y)),1.5)))) * (Foam::log(1 + U/(1-U))));
        uw.x() = (Z / (2*Rc)) * (y-2*Z*y*y0*F1);
        uw.y() = (Foam::sqr(Z)/Rc) * (y0*F1);
    }
    else
    {   
        scalar F2 = (1 / (2*Z*y0)) * (Foam::pow(((1+Foam::sqr(y0))/(1+Foam::sqr(y))),0.25)) * ((1/(U-1)) - (1/(24*Z)*((9*Foam::sqr(y0)+2)/(Foam::pow((1+Foam::sqr(y0)),1.5)) + (3*Foam::sqr(y)-2)/(Foam::pow((1+Foam::sqr(y)),1.5)))) * (Foam::log(1 + 1/(U-1))));
        uw.x() = -(Foam::sqr(Z) / (Rc)) * (y*y0*F2);
        uw.y() =  (Z / (2*Rc)) * (1+2*Z*y0*F2);
    }
    return uw;
}
 
void openPropVLM2D::update()
{
    // Update the time step size.
    dt = runTime_.deltaT().value();

    // Update the current simulation time.
    time = runTime_.timeName();

    // Apply the rotor speed controller and compute rotor speed
    // at new time step.  Also, compute the change in rotor
    // azimuth at the new time step, then rotate blades.
    computeRotSpeed();
    rotateBlades();

    // Apply the yaw controller and compute the change in yaw
    // angle at the new time step, then rotate nacelle.
    computeNacYaw();
    yawNacelle();

    // Find out which processors control which actuator points.
    findControlProcNo();

    // Find the wind vector at each actuator point.
    computeWindVectors();
    iterateInducedVelocities();

    computeVSTARmag();

    // Compute the lift and drag at each actuator point at the 
    // new time step.
    computeBladeForce();

    // Translate the actuator point forces in to flow field body
    // forces.
    computeBodyForce();

}

     
void openPropVLM2D::printDebug()
{
    Info << "Print Debugging Information" << endl;
    Info << "propellerType = " << propellerType << endl;
    Info << "baseLocation = " << baseLocation << endl;
    Info << "numBladePoints = " << numBladePoints << endl;
    Info << "pointDistType = " << pointDistType << endl;
    Info << "epsilon = " << epsilon << endl;
    Info << "smearRadius = " << smearRadius << endl;
    Info << "sphereRadiusScalar = " << sphereRadiusScalar << endl;
    Info << "azimuth = " << azimuth << endl;
    Info << "rotSpeed = " << rotSpeed << endl;
    Info << "nacYaw = " << nacYaw << endl << endl << endl;
    
    Info << "numPropellersDistinct = " << numPropellersDistinct << endl;
    Info << "propellerTypeDistinct = " << propellerTypeDistinct << endl;
    Info << "propellerTypeID = " << propellerTypeID << endl << endl << endl;;

    Info << "NumBl = " << NumBl << endl;
    Info << "TipRad = " << TipRad << endl;
    Info << "HubRad = " << HubRad << endl;
    Info << "OverHang = " << OverHang << endl;
    Info << "ShftTilt = " << ShftTilt << endl;
    Info << "Rake = " << Rake << endl;
    Info << "YawRate = " << YawRate << endl;
    Info << "SpeedControllerType = " << SpeedControllerType << endl;
    Info << "YawControllerType = " << YawControllerType << endl;
    Info << "sectionPoints = " << sectionPoints << endl;
   // Info << "startPt = " << startPt << endl;
    //Info << "endPt = " << endPt << endl;
    Info << "BladeData = " << BladeData << endl;
    Info << "BladeStation = " << BladeStation << endl;
    Info << "BladeChord = " << BladeChord << endl;
    Info << "camberSlope = " << camberSlope << endl;

    Info << "sphereCells = " << sphereCells << endl << endl << endl;

    Info << "db = " << db << endl;
    Info << "bladePoints = " << bladePoints << endl;
    Info << "bladeRadius = " << bladeRadius << endl;
    Info << "boss = " << boss << endl;
    Info << "rotorApex = " << rotorApex << endl;
    Info << "uvShaft = " << uvShaft << endl;
    Info << "uvShaftDir = " << uvShaftDir << endl;
    Info << "uvZ = " << uvZ << endl;
    Info << "deltaNacYaw = " << deltaNacYaw << endl;
    Info << "deltaAzimuth = " << deltaAzimuth << endl;

    Info << "bladeForce = " << bladeForce << endl;
    Info << "windVectors = " << windVectors << endl;
    Info << "bladeAlignedVectors = " << bladeAlignedVectors << endl;
}


volVectorField& openPropVLM2D::force()
{
    // Return the body force field to the solver
    return bodyForce;
}
// Matrix member function
//Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::FillMatrix(Foam::scalarSquareMatrix& M, label n)
Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::FillMatrix(Foam::scalarRectangularMatrix& M, label n)
{
label i,j;
scalarRectangularMatrix R(1,n*n);
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
R[0][n*i+j]=M[i][j];
}
}
return R;
}

Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::InverseMatrix(scalarRectangularMatrix& mAin,label n)
{
scalarRectangularMatrix mAout(1,n*n);
label i,j,k;
scalar sum;
for(i=0;i<n*n;i++){
mAout[0][i]=mAin[0][i];
}

for(i=0;i<n;i++)
{
j=i*n+i;
if((mAout[0][j]<1e-12)&&(mAout[0][j]>-1e-12)){mAout[0][j]=1e-12;}
}

for(i=1;i<n;i++)
{
mAout[0][i]/=(mAout[0][0]);
}

for (i=1; i < n; i++)
{
for (j=i; j < n; j++)
{
sum = 0.0;
for (k = 0; k < i; k++)
{
label index1 =j*n+k;
label index2 =k*n+i;
sum += mAout[0][index1] * mAout[0][index2];
}
label index3 =j*n+i;
mAout[0][index3] -= sum;
}
if (i == n-1) continue;
for (j=i+1.0; j < n; j++)
{
sum = 0.0;
for (k = 0; k < i; k++)
{
label index4=i*n+k;
label index5=k*n+j;
sum += mAout[0][index4]*mAout[0][index5];
}
label index6=i*n+j;
label index7=i*n+i;
mAout[0][index6] = (mAout[0][index6]-sum)/ mAout[0][index7];
}

}

for (i = 0; i < n; i++ )
{
for ( j = i; j < n; j++ )
{
scalar x = 1.0;
if ( i != j )
{
x = 0.0;
for ( k = i; k < j; k++ )
{
label index8=j*n+k;
label index9=k*n+i;
x -= mAout[0][index8]*mAout[0][index9];
}
}
label index10=j*n+i;
label index11=j*n+j;
mAout[0][index10] = x / mAout[0][index11];
}
}

for ( i = 0; i <n; i++ )
{
for ( j = i; j < n; j++ )
{
if ( i == j ) continue;
sum = 0.0;
for ( k = i; k < j; k++ )
{
label index12=k*n+j;
label index13=i*n+k;
sum += mAout[0][index12]*( (i==k) ? 1.0 : mAout[0][index13] );
}
label index14=i*n+j;
mAout[0][index14] = -sum;
}
}

for ( i = 0; i < n; i++ )
{
for ( j = 0; j <n; j++ )
{
sum = 0.0;
for ( k = ((i>j)?i:j); k < n; k++ )
{
label index15=j*n+k;
label index16=k*n+i;
sum += ((j==k)?1.0:mAout[0][index15])*mAout[0][index16];
}
label index17=j*n+i;
mAout[0][index17] = sum;
}
}
return mAout;
}

Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::GetMatrix(Foam::scalarRectangularMatrix& M, label n)
//Foam::scalarSquareMatrix Foam::propellerModels::openPropVLM2D::GetMatrix(Foam::scalarRectangularMatrix& M, label n)
{
label i,j;
scalarRectangularMatrix s(n,n,0);
//scalarSquareMatrix s(n,0);
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
s[i][j]=M[0][n*i+j];
}
}
return s;
}

Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::MultiplyMatrixField
(Foam::scalarRectangularMatrix& S, scalarRectangularMatrix& f, label n, label m)
//(Foam::scalarSquareMatrix& S, scalarRectangularMatrix& f, label n, label m)
{
label i,j,k;
scalarRectangularMatrix Mf(n,m,0);
for(i=0;i<n;i++)
{
for(j=0;j<m;j++)
{
Mf[i][j]=0;
for(k=0;k<n;k++)
{
Mf[i][j]+=S[i][k]*f[k][j];
}
}
}
return Mf;
}

Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::sum_rows
(Foam::scalarRectangularMatrix& S, label n, label m)
{
label i,j;
scalarRectangularMatrix Sr(1,m,0);
for (i=0;i<n;i++)
{
for (j=0;j<m;j++)
{
Sr[i][j] += S[i][j];
}
}
return Sr;
}

Foam::scalarRectangularMatrix Foam::propellerModels::openPropVLM2D::MultiplyField
(Foam::scalarRectangularMatrix& M, scalarField& f, label m)
{
label i,j,k;
scalarRectangularMatrix Mf(1,m,0.0);
for (i=0;i<1;i++)
{
for (j=0;j<m;j++)
{
Mf[i][j] = 0;
for (k=0;k<m;k++)
{
Mf[i][j] += M[i][k]*f[k];
}
}
}
return Mf;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace propellerModels
} // End namespace Foam

// ************************************************************************* //

