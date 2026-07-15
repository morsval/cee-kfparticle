
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
#ifndef CEET0HIT_H
#define CEET0HIT_H 1


#include "FairHit.h"

#include "TObject.h"
#include "TVector3.h"

      /**
       ** @brief Constructor with detailled assignment.
       ** @param[in] detID   Rpc Id. 
       ** @param[in] pos     hit position. 
       ** @param[in] dpos    hit position error.
       ** @param[in] t0[2]  Time Over Threshold [ps].
       ** @param[in] t0[0]  hit time 
       **/
class CeeT0Hit : public FairHit
{

  public:

    /** Default constructor **/
    CeeT0Hit();

    CeeT0Hit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index, const Double_t t0[1]);


    /** Destructor **/
    virtual ~CeeT0Hit();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    Double_t GetTime() { return m_t0[0]; }
    const Double_t* GetT0() const { return m_t0; }
    Double_t GetTriTime() {return m_t0[0];}

  public:

  private:
    Double_t  m_t0[1];

    /** Copy constructor **/
    CeeT0Hit(const CeeT0Hit& hit);
    CeeT0Hit operator=(const CeeT0Hit& hit);

    ClassDef(CeeT0Hit,1)

};

#endif
