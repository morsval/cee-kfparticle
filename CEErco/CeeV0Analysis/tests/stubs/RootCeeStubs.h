#ifndef CEE_V0_ANALYSIS_ROOT_CEE_STUBS_H
#define CEE_V0_ANALYSIS_ROOT_CEE_STUBS_H

// Compile-shape stubs only. They let CI without ROOT/FairRoot parse the full
// analysis translation unit; they are never linked into the real executable.

#include <vector>

typedef long long Long64_t;

class TObject
{
 public:
  virtual ~TObject() {}
};

class TBranch : public TObject {};
class TObjArray : public TObject {};
class TList : public TObject {};

class TAxis
{
 public:
  void SetBinLabel(int, const char*) {}
};

class TIter
{
 public:
  explicit TIter(TList*) {}
  TObject* operator()() { return 0; }
};

class TSystemFile : public TObject
{
 public:
  const char* GetName() const { return ""; }
  bool IsDirectory() const { return false; }
};

class TSystemDirectory
{
 public:
  TSystemDirectory(const char*, const char*) {}
  TList* GetListOfFiles() { return 0; }
};

class TClonesArray : public TObject
{
 public:
  explicit TClonesArray(const char*) {}
  int GetEntriesFast() const { return 0; }
  int GetEntries() const { return 0; }
  TObject* At(int) const { return 0; }
};

class TChain : public TObject
{
 public:
  explicit TChain(const char*) {}
  int Add(const char*) { return 1; }
  Long64_t GetEntries() const { return 0; }
  Long64_t GetEntry(Long64_t) { return 1; }
  TBranch* GetBranch(const char*) const { return 0; }
  void SetBranchStatus(const char*, int) {}
  int SetBranchAddress(const char*, TClonesArray**) { return 0; }
};

class TDirectory { public: class TContext { public: ~TContext() {} }; };

class TFile : public TObject
{
 public:
  static TFile* Open(const char*, const char*) { return 0; }
  bool IsZombie() const { return false; }
  void cd() {}
  int Write() { return 1; }
  TObject* Get(const char*) { return 0; }
  enum { kWriteError = 1 };
  bool TestBit(unsigned) const { return false; }
  void Close() {}
};

class TTree : public TObject
{
 public:
  TTree(const char*, const char*) {}
  template <class T>
  void Branch(const char*, T*, const char*) {}
  int Fill() { return 1; }
  void SetDirectory(TFile*) {}
  TBranch* GetBranch(const char*) { return 0; }
};

class TH1D : public TObject
{
 public:
  TH1D(const char*, const char*, int, double, double) {}
  void Fill(double) {}
  TAxis* GetXaxis() { return &fAxis; }
 private:
  TAxis fAxis;
};

class TNamed : public TObject
{
 public:
  TNamed(const char*, const char*) {}
  void Write() {}
};

class TVector3
{
 public:
  TVector3() : fX(0.), fY(0.), fZ(0.) {}
  double X() const { return fX; }
  double Y() const { return fY; }
  double Z() const { return fZ; }
 private:
  double fX;
  double fY;
  double fZ;
};

class TMatrixDSym
{
 public:
  TMatrixDSym() : fRows(0) {}
  int GetNrows() const { return fRows; }
  int GetNcols() const { return fRows; }
  void ResizeTo(int rows, int) { fRows = rows; }
  double operator()(int row, int column) const { return fValues[row][column]; }
  double& operator()(int row, int column) { return fValues[row][column]; }
 private:
  int fRows;
  double fValues[6][6];
};

class CeeTpcTrack : public TObject
{
 public:
  bool HasVertexSeed() const { return false; }
  const TVector3& GetVertexSeedPos() const { return fVector; }
  const TVector3& GetVertexSeedMom() const { return fVector; }
  void GetVertexSeedCov(TMatrixDSym&) const {}
  bool GetIsFitConverged() const { return true; }
  double GetX() const { return 0.; }
  double GetY() const { return 0.; }
  double GetZ() const { return 0.; }
  double GetPx() const { return 0.; }
  double GetPy() const { return 0.; }
  double GetPz() const { return 1.; }
  double GetFitCharge() const { return 0.; }
  double GetFitChi2() const { return 0.; }
  double GetFitChi2Ndf() const { return 0.; }
  double GetFitNdf() const { return 1.; }
  int GetnHitsFit() const { return 0; }
  double GetdEdx() const { return 0.; }
  int GetMCTrackID() const { return -1; }
 private:
  TVector3 fVector;
};

class CeeTpcVertex : public TObject
{
 public:
  const std::vector<double>& GetX() const { return fDouble; }
  const std::vector<double>& GetY() const { return fDouble; }
  const std::vector<double>& GetZ() const { return fDouble; }
  const std::vector<double>& GetXErr() const { return fDouble; }
  const std::vector<double>& GetYErr() const { return fDouble; }
  const std::vector<double>& GetZErr() const { return fDouble; }
  const std::vector<int>& GetRank() const { return fInt; }
 private:
  std::vector<double> fDouble;
  std::vector<int> fInt;
};

class CeeT0Hit : public TObject
{
 public:
  double GetTime() { return 0.; }
};

class CeeiTOFHit : public TObject
{
 public:
  int GetDetectorID() const { return 0; }
  double GetX() const { return 0.; }
  double GetY() const { return 0.; }
  double GetZ() const { return 0.; }
  double GetTime() const { return 0.; }
};

class CeeEtofHit : public TObject
{
 public:
  double GetX() const { return 0.; }
  double GetY() const { return 0.; }
  double GetZ() const { return 0.; }
  double GetTime() const { return 0.; }
};

class CeeMCTrack : public TObject
{
 public:
  int GetPdgCode() const { return 0; }
  int GetMotherId() const { return -1; }
  double GetPt() const { return 0.; }
  double GetRapidity() const { return 0.; }
};

#endif
