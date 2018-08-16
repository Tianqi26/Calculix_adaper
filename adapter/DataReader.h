/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef DATAREADER_H
#define DATAREADER_H

typedef struct InterfaceData {
	char * facesMeshName;
	char * nodesMeshName;
	char * patchName;
	int numWriteData;
	int numReadData;
	char ** writeDataNames;
	char ** readDataNames;
	double  Calculate_dt;
	double  Calculate_T;
} InterfaceData;

void DataReader_Read(char * dataFilename, char * participantName, char ** CalxfoamDataFilename, InterfaceData ** interfaces, int * numInterfaces);


#endif
