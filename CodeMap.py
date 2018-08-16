##Algorihtm Reference
1. https://github.com/precice/precice/wiki/Adapter-Example
2. https://www5.in.tum.de/pub/Chourdakis2017_Thesis.pdf @ around Page 45
3. https://mediatum.ub.tum.de/doc/1320661/document.pdf

Main=ccx2.13.c
 *func=nonlingeo_precice

->nonlingeo_precice.c
  **1 struct SimulationData simulationData = {
      .ialset = ialset,
      .ielmat = ielmat,
      .istartset = istartset,
      .iendset = iendset,
      .kon = kon,
      .ipkon = ipkon,
      .lakon = &lakon,
      .co = co,
      .set = set,
      .nset = *nset,
      .ikboun = ikboun,
      .ikforc = ikforc,
      .ilboun = ilboun,
      .ilforc = ilforc,
      .nboun = *nboun,###
      .nforc = *nforc,###
      .nelemload = nelemload,
      .nload = *nload,###
      .sideload = sideload,
      .mt = mt,
      .nk = *nk,###
      .theta = &theta,
      .dtheta = &dtheta,
      .tper = tper,
      .nmethod = nmethod,
      .xload = xload,
      .xforc = xforc,
      .xboun = xboun,
      .ntmat_ = ntmat_,
      .vold = vold,
      .fn = fn,
      .cocon = cocon,
      .ncocon = ncocon,
      .mi = mi
   };
   **2 Precice_Setup( configFilename, preciceParticipantName, &simulationData );
   **3 while( Precice_IsCouplingOngoing() ){
       ***3.1 Precice_AdjustSolverTimestep( &simulationData );
       ***3.2 Precice_ReadCouplingData( &simulationData );
          ->adapter/PreciceInterface.c(h)
          ****case FORCES:
		  ****// Read and set forces as concentrated loads (Neumann BC)
          ###Important
		  ****precicec_readBlockVectorData( interfaces[i]->forcesDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
		  ****setNodeForces( interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData, interfaces[i]->numNodes, interfaces[i]->xforcIndices, sim->xforc)
       }
    **4 if( Precice_IsWriteCheckpointRequired() )
        *** 4.1 Precice_WriteIterationCheckpoint( &simulationData, vini );
            ->adapter/PreciceInterface.c(h)
            **** // Save solution vector v
	        **** memcpy( sim->coupling_init_v, v, sizeof( double ) * sim->mt * sim->nk );
        *** 4.2 Precice_FulfilledWriteCheckpoint();
    **5 // Adapter: Perform coupling related actions, only if solver iterations converged (icutb == 0) */
        if( icutb == 0 ){
           ###Important
           *** 5.1 Precice_WriteCouplingData( &simulationData );
               ->adapter/PreciceInterface.c(h)
               ****case DISPLACEMENTS:
			   ****getNodeDisplacements( interfaces[i]->nodeIDs, interfaces[i]->numNodes, sim->vold, sim->mt, interfaces[i]->nodeVectorData );
			   ****precicec_writeBlockVectorData( interfaces[i]->displacementsDataID, interfaces[i]->numNodes, interfaces[i]->preciceNodeIDs, interfaces[i]->nodeVectorData );
           *** 5.2 Precice_Advance( &simulationData );
           *** 5.3 if( Precice_IsReadCheckpointRequired() ){
                        if( *nmethod == 4 )
                       {
                        Precice_ReadIterationCheckpoint( &simulationData, vold );
                        icutb++;
                       }
                       Precice_FulfilledReadCheckpoint();
                    }//if IsRead

        }//if icutb
    **6  Precice_FreeData( &simulationData );


Simplifed Structure:
1. struct SimulationData simulationData
2. CFD Computation
3. Precice_ReadCouplingData( &simulationData )
4. Calculix Computation
5. Precice_WriteCouplingData( &simulationData )
6. Next Timestep and Repeate 2-5

Notes:
1. All func name with Precice_** is comes from the calculix adapter(e.g. adapter/PreciceInterface.h)
2. All func name with precicec_** is comes from the Precice Library(#include "precice/adapters/c/SolverInterfaceC.h")

