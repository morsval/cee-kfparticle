/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *              GNU Lesser General Public Licence (LGPL) version 3,             *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#include "CeeTpcVertex.h"
ClassImp(CeeTpcVertex)

CeeTpcVertex::CeeTpcVertex()
{
  // Default constructor - vectors are automatically initialized as empty
}

CeeTpcVertex::~CeeTpcVertex()
{
  // Destructor - vectors are automatically cleaned up
  Clear();
}

void CeeTpcVertex::AddTrackDCA(int vertexIdx, int trackId, double dcaXY, double dcaZ, double dca3D)
{
  fTrackIds.push_back(trackId);
  fTrackVertexIdx.push_back(vertexIdx);
  fTrackDcaXY.push_back(dcaXY);
  fTrackDcaZ.push_back(dcaZ);
  fTrackDca3D.push_back(dca3D);
}



