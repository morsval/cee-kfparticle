
/********************************************************************************
 *    Copyright (C) 2019 Quark Matter Research Center, IMP,CAS			*
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
/** @author Yingjie Zhou <zhouyingjie@mail.ustc.edu.cn>
 ** @date 2024.01.08
 **/
#include "CeeT0Hit.h"


#include <iostream>
using std::cout;
using std::endl;


// -----   Default constructor   -------------------------------------------
CeeT0Hit::CeeT0Hit()
  : FairHit()
{
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
CeeT0Hit::CeeT0Hit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index, const Double_t t0[1])
  : FairHit(detID, pos, dpos, index)
{
  memcpy(m_t0,t0,sizeof(m_t0));
  this->SetTimeStamp(t0[0]);
  this->SetTimeStampError(0.04);
}


// -----   Destructor   ----------------------------------------------------
CeeT0Hit::~CeeT0Hit() { }
// -------------------------------------------------------------------------


// -----   Public method Print   -------------------------------------------
void CeeT0Hit::Print(const Option_t* /*opt*/) const
{
    std::cout<<"DetID : "<<fDetectorID<<"position : ( "<<fX<<" , "<<fY<<" , "<<fZ<<" )  Time : "<<fTimeStamp<<std::endl;
}
// -------------------------------------------------------------------------

ClassImp(CeeT0Hit)

