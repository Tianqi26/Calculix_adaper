/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using Calxfoam         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include <stdlib.h>
#include "CalxfoamInterface.h"
#include "DataReader.h"
#include "unistd.h"


void Calxfoam_Setup( char * dataFilename, char * participantName, CalcSimulationData * Calcsim )
{
	printf( "Setting up Calxfoam participant %s, using Calxfoam file: %s\n", participantName, dataFilename );
	fflush( stdout );

	int i;
	char * CalxfoamDataFilename;
	int number;
	InterfaceData * interfaces;

	// Read the YAML config file
	DataReader_Read( dataFilename, participantName, &CalxfoamDataFilename, &interfaces, &Calcsim->numCalxfoamInterfaces );
	
	fflush( stdout );

	// Create the solver interface and configure it - Alex: Calculix is always a serial participant (MPI size 1, rank 0)
	//precicec_createSolverInterface( participantName, preciceConfigFilename, 0, 1 );

	// Create interfaces as specified in the config file
	Calcsim->CalxfoamInterfaces = (struct CalxfoamInterface**) malloc( Calcsim->numCalxfoamInterfaces * sizeof( CalxfoamInterface* ) );
    
	for( i = 0 ; i < Calcsim->numCalxfoamInterfaces ; i++ )
	{
		Calcsim->CalxfoamInterfaces[i] = malloc( sizeof( CalxfoamInterface ) );

		CalxfoamInterface_Create( Calcsim->CalxfoamInterfaces[i], Calcsim, &interfaces[i] );
	}
	
	
	// Initialize variables needed for the coupling
	NNEW( Calcsim->coupling_init_v, double, Calcsim->mt * Calcsim->nk );
    
	
	// Initialize Calxfoam time

	Calxfoam_IntializeTime(Calcsim, &interfaces[0]);

	//InterfaceData *iface0=&interfaces[0];

	//Calcsim->Calxfoam_dt = iface0->Calculate_dt; //!here read data and pause
	
    
    
	// Initialize coupling data
	Calxfoam_InitializeData( Calcsim );
	



}



void Calxfoam_InitializeData( CalcSimulationData * Calcsim )
{
	printf( "Initializing coupling data\n" );
	fflush( stdout );

	Calxfoam_WriteCouplingData( Calcsim );
	
	//Calxfoam_ReadCouplingData( Calcsim );
}

void Calxfoam_IntializeTime(CalcSimulationData * Calcsim, InterfaceData * Data)
{
	Calcsim->Calxfoam_dt = Data->Calculate_dt;
}



void Calxfoam_AdjustSolverTimestep( CalcSimulationData * Calcsim )
{
	if( isSteadyStateSimulation( Calcsim->nmethod ) )
	{
		printf( "Adjusting time step for steady-state step\n" );
		fflush( stdout );

		// For steady-state simulations, we will always compute the converged steady-state solution in one coupling step
		*Calcsim->theta = 0;
		*Calcsim->tper = 1;
		*Calcsim->dtheta = 1;

		// Set the solver time step to be the same as the coupling time step
		Calcsim->solver_dt = Calcsim->Calxfoam_dt;
	}
	else
	{
		printf( "Adjusting time step for transient step\n" );
		printf( "Calxfoam_dt dtheta = %f, dtheta = %f, solver_dt = %f\n", Calcsim->Calxfoam_dt / *Calcsim->tper, *Calcsim->dtheta, fmin( Calcsim->Calxfoam_dt, *Calcsim->dtheta * *Calcsim->tper ) );
		fflush( stdout );

		// Compute the normalized time step used by CalculiX
		*Calcsim->dtheta = fmin( Calcsim->Calxfoam_dt / *Calcsim->tper, *Calcsim->dtheta );

		// Compute the non-normalized time step used by preCICE
		Calcsim->solver_dt = ( *Calcsim->dtheta ) * ( *Calcsim->tper );
	}
}

/*

void Precice_Advance( SimulationData * sim )
{
	printf( "Adapter calling advance()...\n" );
	fflush( stdout );

	sim->precice_dt = precicec_advance( sim->solver_dt );
}

bool Precice_IsCouplingOngoing()
{
	return precicec_isCouplingOngoing();
}

bool Precice_IsReadCheckpointRequired()
{
	return precicec_isActionRequired( "read-iteration-checkpoint" );
}

bool Precice_IsWriteCheckpointRequired()
{
	return precicec_isActionRequired( "write-iteration-checkpoint" );
}

void Precice_FulfilledReadCheckpoint()
{
	precicec_fulfilledAction( "read-iteration-checkpoint" );
}

void Precice_FulfilledWriteCheckpoint()
{
	precicec_fulfilledAction( "write-iteration-checkpoint" );
}

void Precice_ReadIterationCheckpoint( SimulationData * sim, double * v )
{

	printf( "Adapter reading checkpoint...\n" );
	fflush( stdout );

	// Reload time
	*( sim->theta ) = sim->coupling_init_theta;

	// Reload step size
	*( sim->dtheta ) = sim->coupling_init_dtheta;

	// Reload solution vector v
	memcpy( v, sim->coupling_init_v, sizeof( double ) * sim->mt * sim->nk );
}

void Precice_WriteIterationCheckpoint( SimulationData * sim, double * v )
{

	printf( "Adapter writing checkpoint...\n" );
	fflush( stdout );

	// Save time
	sim->coupling_init_theta = *( sim->theta );

	// Save step size
	sim->coupling_init_dtheta = *( sim->dtheta );

	// Save solution vector v
	memcpy( sim->coupling_init_v, v, sizeof( double ) * sim->mt * sim->nk );
}
*/


void Calxfoam_ReadCouplingData( CalcSimulationData * Calcsim )
{

	printf( "Adapter reading coupling data...\n" );
	fflush( stdout );

	CalxfoamInterface ** interfaces = Calcsim->CalxfoamInterfaces;
	int numInterfaces = Calcsim->numCalxfoamInterfaces;
	int i;
	//int j;
	//FILE *fileWriter;

    for( i = 0 ; i < numInterfaces ; i++ )
		{
			interfaces[i]->nodeForces = malloc( interfaces[i]->numNodes * 3 * sizeof( double ) );
			Calxfoam_readPresentForce( interfaces[i],Calcsim);
			// !important: the function should be used latter
			setNodeForces( interfaces[i]->CalxfoamNodeIDs, interfaces[i]->nodeForces, interfaces[i]->numNodes, interfaces[i]->xforcIndices, Calcsim->xforc);
		    
		}
	/*
	fileWriter=fopen("xforceIndex.py", "w");
	fprintf(fileWriter, "ID\txforceIndex\t\n"); 
	for ( i = 0 ; i < numInterfaces ; i++ ) 
	{
		for ( j=0 ; j<interfaces[i]->numNodes; j++)
		{
			fprintf(fileWriter, "%d[%d\t%d\t%d]\n", j,
			interfaces[i]->xforcIndices[3 * j],interfaces[i]->xforcIndices[3 * j+1],interfaces[i]->xforcIndices[3 * j+2]);
   
		}
		
	}
	fclose(fileWriter);
	*/

}

char *Str_cat_Num(char str[],int num,char str_end[]) {
	/*Create a C_sytle char by concatenating char and numbers
	Usage:     
		filename = Str_cat_Num("MiuPermGB_Paraview_Timesteps", (int)(t + 1), ".vtk");
	Str_cat_Num("Species",id) return Species[id]

	Programmer: Bin Wang (binwang.0213@gmail.com),Feng Yin (yin.feng@louisiana.edu)
	Creation:   May, 2017
	*/
	static char str1[256];
	static char str2[256];

	strcpy(str2, str);
	sprintf(str1, "%d", num);
	strcat(str2, str1);
	if (str_end!=NULL) strcat(str2, str_end);

	return str2;
}



void Calxfoam_readPresentForce(CalxfoamInterface * interface, CalcSimulationData * Calcsim)
{
	int i;
	FILE *fileReader;
	FILE *fileWriter;
	char *filename; 
	

	filename = Str_cat_Num("CalxNodeForces", *Calcsim->iinc, ".py");
	printf("the read file name: %s\n",filename);
	if (Calxfoam_isReadDataAvailable(filename))
	{

		fileReader=fopen(filename, "r");
		printf("Read force from file"); 
		for ( i = 0 ; i < interface->numNodes ; i++ ) 
		{
			int nodeID = interface->CalxfoamNodeIDs[i];
			// x-component
			fscanf(fileReader, "%lf", &interface->nodeForces[3 * nodeID]);
			// y-component
			fscanf(fileReader, "%lf", &interface->nodeForces[3 * nodeID + 1]);
			// z-component
			fscanf(fileReader, "%lf", &interface->nodeForces[3 * nodeID + 2]);
		}
		fclose(fileReader);
		//printf("-------------Press enter to continue ...----------"); 
		fileWriter=fopen("CalxWriteReadNodeForces.py", "w");
		fprintf(fileWriter, "The data is from reader\n"); 
		for ( i = 0 ; i < interface->numNodes ; i++ ) 
		{
			fprintf(fileWriter, "[%lf\t%lf\t%lf]\n", 
			interface->nodeForces[3 * i],interface->nodeForces[3 * i+1],interface->nodeForces[3 * i+2]);
		}
		fclose(fileWriter);

	}
	

}

int Calxfoam_isReadDataAvailable(char* filename)
{
	int i = 0;
	if( access( filename, F_OK ) != -1 ) {
		i = 1;

    // file exists
	} else {
    printf("file doesn't exist");
	}

	return i;

}





void Calxfoam_WriteCouplingData( CalcSimulationData * Calcsim )
{

	printf( "Adapter writing coupling data...\n" );
	fflush( stdout );

	CalxfoamInterface ** interfaces = Calcsim->CalxfoamInterfaces;
	int numInterfaces = Calcsim->numCalxfoamInterfaces;
	int i;
	int iset;

	//!initialize value
	//Calcsim->solver_dt == 0.1;
	
	//interfaces[i]->writeData == DISPLACEMENTS;
    for( i = 0 ; i < numInterfaces ; i++ )
	{
		//! Open Memory for nodeDisplacements, remember to free after it is no longer used
		interfaces[i]->nodeDisplacements = malloc( interfaces[i]->numNodes * 3 * sizeof( double ) );
		getNodeDisplacements( interfaces[i]->nodeIDs, interfaces[i]->numNodes, Calcsim->vold, Calcsim->mt, interfaces[i]->nodeDisplacements );
        Calxfoam_writePresentCoordinate( interfaces[i], Calcsim);
		//printf("-----------------the mt of calcsim is%d-------------------\n",Calcsim->mt);
		//Calxfoam_writePresentCoordinate( interfaces[i], Calcsim);
	}


}


void Calxfoam_writePresentCoordinate( CalxfoamInterface * interface, CalcSimulationData * Calcsim  )
{
	int i;
	FILE *fileWriter;
	
	char *filename; 
	

	filename = Str_cat_Num("WriteNewCoordinates", *Calcsim->iinc, ".py");
	printf("the read file name: %s\n",filename);
	fileWriter=fopen(filename, "w");
	fprintf(fileWriter, "The data is Coordinates in the present step\n"); 
	fprintf(fileWriter, "The increment %d\t the total time %lf\n",*Calcsim->iinc, *Calcsim->theta**Calcsim->tper); 
	fprintf(fileWriter, "NodeNumber %d\n",interface->numNodes); 
	fprintf(fileWriter,"Index\tx\ty\tz\n");

	for( i = 0 ; i < interface->numNodes ; i++ )
	{
		fprintf(fileWriter, "%d\t%lf\t%lf\t%lf\n", interface->CalxfoamNodeIDs[i],
			    interface->nodeCoordinates[i * 3 + 0]+interface->nodeDisplacements[i * 3 + 0],interface->nodeCoordinates[i * 3 + 1]+interface->nodeDisplacements[i * 3 + 1],interface->nodeCoordinates[i * 3 + 2]+interface->nodeDisplacements[i * 3 + 2]);
	}

	
	fclose(fileWriter);
}



void Calxfoam_FreeData( CalcSimulationData * Calcsim )
{
	int i;

	if( Calcsim->coupling_init_v != NULL ){
		free( Calcsim->coupling_init_v );
	}

	for( i = 0 ; i < Calcsim->numCalxfoamInterfaces ; i++ )
	{
		CalxfoamInterface_FreeData( Calcsim->CalxfoamInterfaces[i] );
		if( Calcsim->CalxfoamInterfaces[i] != NULL ){
			free( Calcsim->CalxfoamInterfaces[i] );
		}
	}

	
}

//CalxfoamInterface_Create( Calcsim->CalxfoamInterfaces[i], Calcsim, &interfaces[i] );
void CalxfoamInterface_Create( CalxfoamInterface * interface, CalcSimulationData * Calcsim, InterfaceData * Data )
{

	// Initialize pointers as NULL
	interface->elementIDs = NULL;
	interface->faceIDs = NULL;
	interface->faceCenterCoordinates = NULL;
	interface->CalxfoamFaceCenterIDs = NULL;
	interface->nodeCoordinates = NULL;
	interface->CalxfoamNodeIDs = NULL;
	interface->triangles = NULL;
	interface->nodeScalarData = NULL;
	interface->nodeVectorData = NULL;
	interface->faceCenterData = NULL;
	interface->xbounIndices = NULL;
	interface->xloadIndices = NULL;
	interface->xforcIndices = NULL;

	// The patch identifies the set used as interface in Calculix
	interface->name = Data->patchName;

	// Nodes mesh
	interface->nodesMeshID = -1;
	interface->nodesMeshName = Data->nodesMeshName;
	CalxfoamInterface_DataNodesMesh( interface, Calcsim );

	// Face centers mesh
	interface->faceCentersMeshID = -1;
	interface->faceCentersMeshName = Data->facesMeshName;
		//Only configure a face center mesh if necesary; i.e. do not configure it for FSI simulations, also do not configure tetra faces if no face center mesh is used (as in FSI simulations)
		if ( interface->faceCentersMeshName != NULL) {
			//PreciceInterface_ConfigureFaceCentersMesh( interface, sim );
		// Triangles of the nodes mesh (needs to be called after the face centers mesh is configured!)
			//PreciceInterface_ConfigureTetraFaces( interface, sim );
		}

	CalxfoamInterface_ConfigureCouplingData( interface, Calcsim, Data );

}

/*

void PreciceInterface_ConfigureFaceCentersMesh( PreciceInterface * interface, SimulationData * sim )
{

	char * faceSetName = toFaceSetName( interface->name );
	interface->faceSetID = getSetID( faceSetName, sim->set, sim->nset );
	interface->numElements = getNumSetElements( interface->faceSetID, sim->istartset, sim->iendset );

	interface->elementIDs = malloc( interface->numElements * sizeof( ITG ) );
	interface->faceIDs = malloc( interface->numElements * sizeof( ITG ) );
	getSurfaceElementsAndFaces( interface->faceSetID, sim->ialset, sim->istartset, sim->iendset, interface->elementIDs, interface->faceIDs );

	interface->faceCenterCoordinates = malloc( interface->numElements * 3 * sizeof( double ) );
	getTetraFaceCenters( interface->elementIDs, interface->faceIDs, interface->numElements, sim->kon, sim->ipkon, sim->co, interface->faceCenterCoordinates );

	interface->faceCentersMeshID = precicec_getMeshID( interface->faceCentersMeshName );
	interface->preciceFaceCenterIDs = malloc( interface->numElements * sizeof( int ) );
	precicec_setMeshVertices( interface->faceCentersMeshID, interface->numElements, interface->faceCenterCoordinates, interface->preciceFaceCenterIDs );

}
*/
//! get the nodes mesh in the interface and set all the datas
void CalxfoamInterface_DataNodesMesh( CalxfoamInterface * interface, CalcSimulationData * Calcsim )
{

	char * nodeSetName = toNodeSetName( interface->name );
	//char * nodeSetName = "Nsurface";
	int i;
	FILE *fileWriter;

	interface->nodeSetID = getSetID( nodeSetName, Calcsim->set, Calcsim->nset );
	interface->numNodes = getNumSetElements( interface->nodeSetID, Calcsim->istartset, Calcsim->iendset );
	interface->nodeIDs = &Calcsim->ialset[Calcsim->istartset[interface->nodeSetID] - 1]; //Lucia: make a copy

	interface->nodeCoordinates = malloc( interface->numNodes * 3 * sizeof( double ) );
	getNodeCoordinates( interface->nodeIDs, interface->numNodes, Calcsim->co, Calcsim->vold, Calcsim->mt, interface->nodeCoordinates );

	if( interface->nodesMeshName != NULL )
	{
		interface->nodesMeshID = 1;
		interface->CalxfoamNodeIDs = malloc( interface->numNodes * sizeof( int ) );
		Calxfoam_setMeshVertices( interface->nodesMeshID, interface->numNodes, interface->nodeCoordinates, interface->CalxfoamNodeIDs );
	}

	fileWriter=fopen("WriteCoordinates.py", "w");
	fprintf(fileWriter, "The data is from reader\n"); 
	for ( i = 0 ; i < interface->numNodes ; i++ ) 
	{
		fprintf(fileWriter, "[%lf\t%lf\t%lf]\n", 
			    interface->nodeCoordinates[i * 3 + 0],interface->nodeCoordinates[i * 3 + 1],interface->nodeCoordinates[i * 3 + 2]);
	}
	fclose(fileWriter);

}

void Calxfoam_setMeshVertices( int nodesId, int nodesNum, double * CalxfoamnodeCoordinates, int * NodesIDInCalxfoam )
{
	int i;
	int k;
	int j;
	k=0;
	for ( i=0; i<nodesId; i++)
	{
		for (j=0; j<nodesNum; j++)
		{
			NodesIDInCalxfoam[k]=k;
			k++;

		}
	}

}

void CalxfoamInterface_EnsureValidNodesMeshID( CalxfoamInterface * interface )
{
	if( interface->nodesMeshID < 0 )
	{
		printf( "Nodes mesh not provided in YAML config file\n" );
		fflush( stdout );
		exit( EXIT_FAILURE );
	}
}

/*

void PreciceInterface_ConfigureTetraFaces( PreciceInterface * interface, SimulationData * sim )
{
	int i;

	if( interface->nodesMeshName != NULL )
	{
		interface->triangles = malloc( interface->numElements * 3 * sizeof( ITG ) );
		getTetraFaceNodes( interface->elementIDs, interface->faceIDs,  interface->nodeIDs, interface->numElements, interface->numNodes, sim->kon, sim->ipkon, interface->triangles );

		for( i = 0 ; i < interface->numElements ; i++ )
		{
			precicec_setMeshTriangleWithEdges( interface->nodesMeshID, interface->triangles[3*i], interface->triangles[3*i+1], interface->triangles[3*i+2] );
		}
	}
}
*/

//!get xforcIndices

void CalxfoamInterface_ConfigureCouplingData( CalxfoamInterface * interface, CalcSimulationData * Calcsim, InterfaceData * data )
{

	interface->nodeScalarData = malloc( interface->numNodes * sizeof( double ) );
	interface->nodeVectorData = malloc( interface->numNodes * 3 * sizeof( double ) );
	interface->faceCenterData = malloc( interface->numElements * sizeof( double ) );

	int i;

	for( i = 0 ; i < data->numReadData ; i++ )
	{
		
		if( strcmp1( data->readDataNames[i], "Forces" + i ) == 0 )
		{

			CalxfoamInterface_EnsureValidNodesMeshID( interface );
			interface->readData = FORCES;
			interface->xforcIndices = malloc( interface->numNodes * 3 * sizeof( int ) );
			interface->forcesDataID = interface->nodesMeshID;
			getXforcIndices( interface->nodeIDs, interface->numNodes, Calcsim->nforc, Calcsim->ikforc, Calcsim->ilforc, interface->xforcIndices );
			printf( "Read data '%s' found.\n", data->readDataNames[i] );
		}
		else
		{
			printf( "ERROR: Read data '%s' does not exist!\n", data->readDataNames[i] );
			exit( EXIT_FAILURE );
		}
	}

	for( i = 0 ; i < data->numWriteData ; i++ )
	{
		if( strcmp1( data->writeDataNames[i], "Displacements" + i ) == 0 )
		{
			CalxfoamInterface_EnsureValidNodesMeshID( interface );
			interface->writeData = DISPLACEMENTS;
			interface->displacementsDataID = interface->nodesMeshID;
			printf( "Write data '%s' found.\n", data->writeDataNames[i] );
		}
		else
		{
			printf( "ERROR: Write data '%s' does not exist!\n", data->writeDataNames[i] );
			exit( EXIT_FAILURE );
		}
	}
}


void CalxfoamInterface_FreeData( CalxfoamInterface * calxfoamInterface )
{
	if( calxfoamInterface->elementIDs != NULL ){
		free( calxfoamInterface->elementIDs );
	}

	if( calxfoamInterface->faceIDs != NULL ){
		free( calxfoamInterface->faceIDs );
	}

	if( calxfoamInterface->faceCenterCoordinates != NULL ){
		free( calxfoamInterface->faceCenterCoordinates );
	}

	if( calxfoamInterface->CalxfoamFaceCenterIDs != NULL ){
		free( calxfoamInterface->CalxfoamFaceCenterIDs );
	}

	if( calxfoamInterface->nodeCoordinates != NULL ){
		free( calxfoamInterface->nodeCoordinates );
	}

	if( calxfoamInterface->CalxfoamNodeIDs != NULL ){
		free( calxfoamInterface->CalxfoamNodeIDs );
	}

	if( calxfoamInterface->triangles != NULL ){
		free( calxfoamInterface->triangles );
	}

	if( calxfoamInterface->nodeScalarData != NULL ){
		free( calxfoamInterface->nodeScalarData );
	}

	if( calxfoamInterface->nodeVectorData != NULL ){
		free( calxfoamInterface->nodeVectorData );
	}

	if( calxfoamInterface->faceCenterData != NULL ){
		free( calxfoamInterface->faceCenterData );
	}

	if( calxfoamInterface->xbounIndices != NULL ){
		free( calxfoamInterface->xbounIndices );
	}

	if( calxfoamInterface->xloadIndices != NULL ){
		free( calxfoamInterface->xloadIndices );
	}

	if ( calxfoamInterface->xforcIndices != NULL ){
		free( calxfoamInterface->xforcIndices );
	}

	if ( calxfoamInterface->nodeForces != NULL ){
		free( calxfoamInterface->nodeForces );
	}

	if ( calxfoamInterface->nodeDisplacements != NULL ){
		free( calxfoamInterface->nodeDisplacements );
	}


}





void checkCalxfoamInterface(){
	printf ("the CalxfoamInterface Used");
}


