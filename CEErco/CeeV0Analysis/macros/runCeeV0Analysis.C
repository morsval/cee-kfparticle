#include "TString.h"
#include "TSystem.h"
#include "TError.h"

namespace
{
TString ShellQuote(const char* text)
{
  TString output(text ? text : "");
  output.ReplaceAll("'", "'\\''");
  return TString("'") + output + "'";
}
}

// Run from CEErco/CeeV0Analysis after `make`.
// For real data with CeeT0Hit, set t0ScaleToNs explicitly (1 if already ns,
// 1e-3 if the persisted value is ps; verify this with hT0Raw/hT0 first).
int runCeeV0Analysis(const char* input,
                     const char* output = "CeeV0Lambda.root",
                     const char* mode = "simulation",
                     const char* cutPreset = "legacy",
                     const char* t0Source = "",
                     double t0ScaleToNs = -1.0,
                     Long64_t events = -1,
                     const char* executable = "./CeeV0Analysis.bin",
                     Bool_t allowIncompleteMotherJacobian = kFALSE)
{
  const TString dataMode(mode ? mode : "");
  const TString requestedT0(t0Source ? t0Source : "");
  if (dataMode == "simulation" &&
      ((requestedT0.Length() > 0 && requestedT0 != "fixed") ||
       t0ScaleToNs > 0.))
  {
    Error("runCeeV0Analysis",
          "simulation uses exact generator t=0 and cannot read/scale CeeT0Hit");
    return 2;
  }
  TString command = ShellQuote(executable) + " --input " + ShellQuote(input) +
    " --output " + ShellQuote(output) + " --mode " + ShellQuote(mode) +
    " --cut-preset " + ShellQuote(cutPreset) +
    TString::Format(" --events %lld", events);
  if (t0Source && TString(t0Source).Length() > 0)
  {
    command += " --t0-source ";
    command += ShellQuote(t0Source);
  }
  // A negative sentinel means "do not pass the option". In particular this
  // preserves the executable's check that real CeeT0Hit data require the
  // analyst to confirm the raw-time unit explicitly.
  if (t0ScaleToNs > 0.)
    command += TString::Format(" --t0-scale-to-ns %.17g", t0ScaleToNs);
  if (allowIncompleteMotherJacobian)
    command += " --allow-incomplete-mother-jacobian";
  return gSystem->Exec(command);
}
