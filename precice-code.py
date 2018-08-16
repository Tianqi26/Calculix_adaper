Precice_Setup()
{
    ConfigReader_Read()
    {
        class:interfaceConfig;
            1.numInterface;
            2.class interfaceConfig: interfaces
            -nodesMeshName; //Calculix_Mesh
            -facesMeshName;
            -patchName;//surface
            -numWriteData;
            -numReadData;
            -WriteDataNames;//DisplacementDeltas0
            -ReadDataNames;//Forces0
    }

    //!creat interface in precice coupling for changing data

    precicec_createSolverInterface(participantName, preciceConfigFileName,0,1)
    {
        class: SolverInterfaceImpl;
            -accessorName: ParticipantName;
            _accessorProcessRank(accessorProcessRank),
            _accessorCommunicatorSize(accessorCommunicatorSize),
            _accessor(),
            _dimensions(0),
            _serverMode(serverMode),
            _clientMode(false),
            _meshIDs(),
            _dataIDs(),
            _exportVTKNeighbors(),
            _m2ns(),
            _participants(),
            _numberAdvanceCalls(0),
            _requestManager(nullptr)
    }

    struct: sim->preciceInterfaces;

    //!creat and set datas for interface 

    
    PreciceInterface_Create()
    {
        Initialize:
        struct: preciceInterfaces
        
        //node mesh;
        PreciceInterface_ConfigureNodesMesh()
        {
            -nodeSetName: NSURFACEN
            -nodeSetID: 3
            -numNodes: 58
            -nodeIDs
            -nodeCoordinates
            -nodesMeshID: precicec_getMeshID( interface->nodesMeshName )
            -preciceNodeIDs: precicec_setMeshVertices( interface->nodesMeshID, interface->numNodes, interface->nodeCoordinates, interface->preciceNodeIDs )
        }
        //face center mesh;
        PreciceInterface_ConfigureFaceCentersMesh(){}
        PreciceInterface_ConfigureTetraFaces(){}

        // creat the interface and link to the changing data interface
        PreciceInterface_ConfigureCouplingData()
        {
            nodeScalarData: allocate
            nodeVectorData: allocate
            faceCenterData: allocate
            //read data
            1.Temperature: 
              -readData=TEMPERATURE; 
              -temperatureDataID =precicec_getDataID( "Temperature", interface->nodesMeshID )
              -xbounIndices=getXbounIndices( interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, TEMPERATURE );
            2.Heat-Flux: 
            3.Sink-Temperature-
            4.Heat-Transfer-Coefficient-
            5.Forces:
              -readData = FORCES
              -forcesDataID = precicec_getDataID( config->readDataNames[i], interface->nodesMeshID );
              -xforcIndices = getXforcIndices( interface->nodeIDs, interface->numNodes, sim->nforc, sim->ikforc, sim->ilforc, interface->xforcIndices );
            6.Displacements
              -readData = DISPLACEMENTS;
              -displacementsDataID = precicec_getDataID( config->readDataNames[i], interface->nodesMeshID );
              -xbounIndices = getXbounIndices( interface->nodeIDs, interface->numNodes, sim->nboun, sim->ikboun, sim->ilboun, interface->xbounIndices, DISPLACEMENTS );
            
            
            //write data
            1.Temperature: 
              -writeData = TEMPERATURE; 
              -temperatureDataID = precicec_getDataID( "Temperature", interface->nodesMeshID );
            2.Heat-Flux: 
            3.Sink-Temperature-
            4.Heat-Transfer-Coefficient-
            5.Forces:
              -writeData = FORCES;
              -forcesDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
            6.Displacements
              -writeData = DISPLACEMENTS;
              -displacementsDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
            7.DisplacementDeltas
              -writeData = DISPLACEMENTDELTAS;
              -displacementDeltasDataID = precicec_getDataID( config->writeDataNames[i], interface->nodesMeshID );
        }

    }

    allocate coupling_init_v;

    precice_dt = precicec_initialize();



    //!initialize data for interface in precice coupling for changing data

    Precice_InitializeData( sim )
    {
        Precice_WriteCouplingData( sim );
        {
            1.Temperature: 
              -nodeScalarData=getNodeTemperatures( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeScalarData ); 
              -precicec_writeBlockScalarData( interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData );
              
            2.Heat-Flux: 
            3.Sink-Temperature-
            4.Heat-Transfer-Coefficient-
            5.Forces:
              -nodeVectorData=getNodeForces( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->fn, sim->mt, interfaces[i]->nodeVectorData );
			  -precicec_writeBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
            6.Displacements
              -nodeVectorData=getNodeDisplacements( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeVectorData );
			  -precicec_writeBlockVectorData( interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
            7.DISPLACEMENTDELTAS:
              -nodeVectorData=getNodeDisplacementDeltas( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->coupling_init_v, sim->mt, interfaces[i]->nodeVectorData );
              -precicec_writeBlockVectorData( interfaces[i]->displacementDeltasDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
            if( precicec_isActionRequired( "write-initial-data" ) )
            {
                precicec_fulfilledAction( "write-initial-data" );
            }
        }
	    precicec_initialize_data();

        //!set data for interface in calculix; get the data from changing data interface in precice


	    Precice_ReadCouplingData( sim )
        {
            if( precicec_isReadDataAvailable() )
            {
            1. TEMPERATURE:
				// Read and set temperature BC
				nodeScalarData=precicec_readBlockScalarData( interfaces[i]->temperatureDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeScalarData );
				setNodeTemperatures( interfaces[i]->nodeScalarData, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun );
				break;
			2. HEAT_FLUX:
				
			3. CONVECTION:

			4. FORCES:
				// Read and set forces as concentrated loads (Neumann BC)
				nodeVectorData=precicec_readBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
				setNodeForces( interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->xforcIndices, sim->xforc);
				break;
			4. DISPLACEMENTS:
				// Read and set displacements as single point constraints (Dirichlet BC)
				nodeVectorData = precicec_readBlockVectorData( interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
				setNodeDisplacements( interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->xbounIndices, sim->xboun );
				break;
			5. DISPLACEMENTDELTAS:
				printf( "DisplacementDeltas cannot be used as read data\n" );
				fflush( stdout );
				exit( EXIT_FAILURE );
            }
        }
    }

}


while (Precice_IsCouplingOngoing)
{
    Precice_AdjustSolverTimestep( &simulationData )
    {
        theta;
        tper;
        dtheta;
        solver_dt;
        precice_dt
    } //!change control time data from this function

    Precice_ReadCouplingData( &simulationData ) //the 2rd time
    {
        if( precicec_isReadDataAvailable() )
            {
			4. FORCES:
				// Read and set forces as concentrated loads (Neumann BC)
				nodeVectorData=precicec_readBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
				setNodeForces( interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->xforcIndices, sim->xforc);
				break;
            }

    }

    if( Precice_IsWriteCheckpointRequired() )
      {
          Precice_WriteIterationCheckpoint( &simulationData, vini );
          {
              // Save time
              sim->coupling_init_theta = *( sim->theta );
              
              // Save step size
              sim->coupling_init_dtheta = *( sim->dtheta );
              
              // Save solution vector v
              sim->coupling_init_v=vini;
          }
          Precice_FulfilledWriteCheckpoint();
      }


      if( icutb == 0 )
      {
          Precice_WriteCouplingData( &simulationData ); //the 2rd time
          Precice_Advance( &simulationData );
          {
              precice_dt = precicec_advance( sim->solver_dt )
          }
          if( Precice_IsReadCheckpointRequired() )
          {
              if( *nmethod == 4 )
              {
                  Precice_ReadIterationCheckpoint( &simulationData, vold )
                  {
                      theta = coupling_init_theta;
                      dtheta = coupling_init_dtheta;

	                  // Reload solution vector v
	                  vold = coupling_init_v
                  }
                icutb++;
              }
              Precice_FulfilledReadCheckpoint();
          }
      }


}