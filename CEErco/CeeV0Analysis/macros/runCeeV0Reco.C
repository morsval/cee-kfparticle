// CEE v0.6.5 Lambda reconstruction launcher (Linux / ROOT).
// Put this file in CEErco/CeeV0Analysis/macros/ and run from CeeV0Analysis/.
// First compile the existing analysis: make -j4
//
// Single ROOT file, first 10000 events:
// root -l -b -q 'macros/runCeeV0Reco.C("/data/reco.root","lambda.root",10000)'
// ROOT file list, all events:
// root -l -b -q 'macros/runCeeV0Reco.C("/data/reco.list","lambda.root",-1)'
// Magnetic tracking, signed field along y in tesla:
// root -l -b -q 'macros/runCeeV0Reco.C("/data/reco.list","lambda_By.root",-1,"tracking",0.5)'
//
// This macro invokes CeeV0Analysis.bin. Reconstruction and ROOT I/O stay in
// the compiled analysis; it is not an independent track/vertex fitter.
// Existing analysis entry: Lambda -> p + pi-, TPC-TPC only, not hypertriton.
// MC truth is only an optional output label, never a reconstruction cut.

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <string>
#include <sys/wait.h>

namespace CeeV0RecoMacro
{
// Quote each argument separately. Spaces, apostrophes and shell metacharacters
// in file paths are literal; do not expand ~ or environment variables here.
inline std::string Quote(const std::string& input)
{
  std::string out = "'";
  for (std::size_t i = 0; i < input.size(); ++i)
    out += input[i] == '\'' ? "'\\''" : std::string(1, input[i]);
  return out + "'";
}

template<class T>
void Option(std::ostringstream& command, const char* key, const T& value)
{
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << std::setprecision(17) << value;
  command << ' ' << key << ' ' << Quote(text.str());
}
}

int runCeeV0Reco(const char* input,
                 const char* output = "lambda.root",
                 long long events = -1,
                 const char* momentum = "tof",
                 double byTesla = 0.,
                 const char* mode = "simulation",
                 const char* executable = "./CeeV0Analysis.bin",
                 bool dryRun = false)
{
  // ===== 修改这里的配置和 cut，不需要改重建核心 =====
  const char* treeName = "ceesim";
  // 有磁默认：沿测量 track 移到与当前重建 PV 的 DCA 点。
  // 同步输运动量及完整 6x6 协方差；不是对子 track 加 PV 约束。
  const char* magneticTrackReference = "pv-dca"; // stored 仅用于旧参考对照
  const double referenceMaximumAbsPathCm = 1000.;
  const int maximumInputFiles = -1;           // -1: 全部文件
  const int debugLevel = 1;                   // 0/1/2
  const int maximumTracksPerEvent = 1000;
  const int maximumCandidatesPerEvent = 100000;

  // track 输入质量。下面沿用宽松入口，数值不是优化后的物理工作点。
  const int minimumTpcHits = 0;               // 例如改为 15
  const double maximumTrackChi2Ndf = 1.e30;   // 例如改为 5.0
  const bool requireTrackConverged = false;
  const bool strictPvErrors = false;          // true: 拒绝缺失/无效的 PV 误差
  const bool disableLegacyTrackFallback = false;
  // 有场 tracking 无论这个开关如何，都必须有有效六维种子及协方差。

  // 零场 track 与 TOF hit 的匹配上限，单位 cm。有场 tracking 不使用这些 hit。
  const double maximumItofMatchDistanceCm = 10.;
  const double maximumEtofMatchDistanceCm = 10.;
  const bool disableDcaRole = false;
  const double lambdaRoleDcaBoundaryCm = 1.2;
  // 零场默认: proton 假设 DCA_to_PV <= 1.2，pion 假设 > 1.2 cm。
  // tracking 模式由分析程序自动关闭这条零场角色规则，使用拟合电荷。

  // Lambda 候选 cut。1.e30 表示数值上放开上限。
  // 这里显式覆盖 legacy 预设。改变这些值才会改变本宏的候选 cut。
  const double minimumMassGeV = 0.7;
  const double maximumMassGeV = 2.0;
  const double maximumDaughterDcaCm = 1.e30;      // 例如 1.0 cm
  const double maximumGeometricChi2Ndf = 1.e30;   // 例如 10.0
  const double maximumTopologicalChi2Ndf = 1.e30; // 几何项+旧 PV 项的总 chi2/NDF
  const double minimumDecayLengthCm = 0.;        // 约束前 |SV-PV|，不是带电弧长
  const double minimumDecayLengthSignificance = 0.; // 约束前长度显著性，先不收紧
  const double minimumPointingCosine = -1.;
  const double minimumDaughterChi2FromPv = 0.;    // 两条子 track 各自的原始 chi2
  // minimumDecayLengthSignificance 对应 --min-ldl，切的是 decayLengthSignificance。
  // 它不是母粒子 PV 更新后的 motherTopoDecayLengthSignificance。
  // maximumTopologicalChi2Ndf 也不是 motherTopoPvChi2/motherTopoPvNdf。
  // 新八维 PV 拟合量请在输出 LambdaCandidates 上再选，见配套 README。

  // 误差参数与时间单位：沿用程序默认值，须按自己的文件/标定核对。
  const double tofRawToNsScale = 1.;
  const double tofOffsetNs = 0.;
  const double tofTimeSigmaNs = 0.10;
  const double tofPathLengthSigmaCm = 0.5;
  const double maximumBeta = 0.9995;
  const double fallbackTrackSigmaXcm = 0.10;
  const double fallbackTrackSigmaYcm = 0.10;
  const double fallbackTrackSigmaZcm = 0.15;
  const double fallbackDirectionSigmaRad = 0.003;
  const double fallbackPvSigmaXcm = 0.10;
  const double fallbackPvSigmaYcm = 0.10;
  const double fallbackPvSigmaZcm = 0.15;

  // 只有 mode="real" 才使用这组 T0 参数。
  // simulation 固定生成器时间原点为 0，不读取 CeeT0Hit，也不加入公共 T0 误差。
  const char* realT0Source = "hit";             // hit 或 fixed
  const double realT0ScaleToNs = -1.;           // 确认单位后填 1(ns) 或 0.001(ps)
  const double realT0OffsetNs = 0.;
  const char* realT0Policy = "median";
  const double realFixedT0Ns = 0.;              // 仅 fixed 使用
  const double realT0SigmaNs = 0.10;            // 占位，按标定修改
  const double realAdditionalTimeSigmaNs = 0.10;// 占位，勿重复计入同一误差

  const bool saveSimulationTruth = true;
  const bool buildProductionState = true;      // 保存独立八维 PV 约束副本
  // 保留默认 TOF 收敛及完整母粒子 Jacobian 要求，不启用诊断性放宽开关。
  // ===== 配置结束 =====

  if (!input || !*input || !output || !*output || !executable || !*executable) {
    std::cerr << "Input, output and executable must be non-empty.\n";
    return 2;
  }
  const std::string mom = momentum ? momentum : "";
  const std::string dataMode = mode ? mode : "";
  if ((mom != "tof" && mom != "tracking") ||
      (dataMode != "simulation" && dataMode != "real") ||
      !std::isfinite(byTesla) || events < -1) {
    std::cerr << "Use tof|tracking, simulation|real, finite By, events >= -1.\n";
    return 2;
  }
  if ((mom == "tof" && std::fabs(byTesla) > 1.e-6) ||
      (mom == "tracking" && std::fabs(byTesla) <= 1.e-6)) {
    std::cerr << "tof requires zero field; tracking requires nonzero By.\n";
    return 2;
  }
  if (dataMode == "real" && std::string(realT0Source) == "hit" &&
      !(realT0ScaleToNs > 0.)) {
    std::cerr << "Set realT0ScaleToNs after checking the persisted T0 unit.\n";
    return 2;
  }

  using CeeV0RecoMacro::Option;
  std::ostringstream command;
  command << CeeV0RecoMacro::Quote(executable);
  Option(command, "--input", input);
  Option(command, "--output", output);
  Option(command, "--tree", treeName);
  Option(command, "--mode", dataMode);
  Option(command, "--momentum", mom);
  Option(command, "--by-tesla", byTesla);
  Option(command, "--magnetic-track-reference", magneticTrackReference);
  Option(command, "--reference-max-path-cm", referenceMaximumAbsPathCm);
  Option(command, "--events", events);
  Option(command, "--files", maximumInputFiles);
  Option(command, "--debug", debugLevel);
  Option(command, "--cut-preset", "legacy");
  Option(command, "--max-tracks", maximumTracksPerEvent);
  Option(command, "--max-candidates", maximumCandidatesPerEvent);
  Option(command, "--min-tpc-hits", minimumTpcHits);
  Option(command, "--max-track-chi2ndf", maximumTrackChi2Ndf);
  Option(command, "--itof-match-cm", maximumItofMatchDistanceCm);
  Option(command, "--etof-match-cm", maximumEtofMatchDistanceCm);
  Option(command, "--lambda-role-dca-cm", lambdaRoleDcaBoundaryCm);
  Option(command, "--mass-min", minimumMassGeV);
  Option(command, "--mass-max", maximumMassGeV);
  Option(command, "--max-daughter-dca", maximumDaughterDcaCm);
  Option(command, "--max-geo-chi2ndf", maximumGeometricChi2Ndf);
  Option(command, "--max-topo-chi2ndf", maximumTopologicalChi2Ndf);
  Option(command, "--min-decay-length", minimumDecayLengthCm);
  Option(command, "--min-ldl", minimumDecayLengthSignificance);
  Option(command, "--min-pointing-cos", minimumPointingCosine);
  Option(command, "--min-daughter-pv-chi2", minimumDaughterChi2FromPv);
  Option(command, "--tof-time-scale-to-ns", tofRawToNsScale);
  Option(command, "--tof-time-offset-ns", tofOffsetNs);
  Option(command, "--tof-time-sigma-ns", tofTimeSigmaNs);
  Option(command, "--tof-path-sigma-cm", tofPathLengthSigmaCm);
  Option(command, "--max-beta", maximumBeta);
  Option(command, "--fallback-track-sigma-x-cm", fallbackTrackSigmaXcm);
  Option(command, "--fallback-track-sigma-y-cm", fallbackTrackSigmaYcm);
  Option(command, "--fallback-track-sigma-z-cm", fallbackTrackSigmaZcm);
  Option(command, "--fallback-direction-sigma-rad", fallbackDirectionSigmaRad);
  Option(command, "--fallback-pv-sigma-x-cm", fallbackPvSigmaXcm);
  Option(command, "--fallback-pv-sigma-y-cm", fallbackPvSigmaYcm);
  Option(command, "--fallback-pv-sigma-z-cm", fallbackPvSigmaZcm);
  if (requireTrackConverged) command << " --require-track-converged";
  if (strictPvErrors) command << " --strict-pv-errors";
  if (disableLegacyTrackFallback) command << " --disable-legacy-track-fallback";
  if (disableDcaRole) command << " --disable-dca-role";
  if (!buildProductionState) command << " --no-production-state";
  if (!saveSimulationTruth) command << " --no-truth";

  if (dataMode == "real") {
    Option(command, "--t0-source", realT0Source);
    Option(command, "--t0-sigma-ns", realT0SigmaNs);
    Option(command, "--additional-decay-time-sigma-ns", realAdditionalTimeSigmaNs);
    if (std::string(realT0Source) == "hit") {
      Option(command, "--t0-scale-to-ns", realT0ScaleToNs);
      Option(command, "--t0-offset-ns", realT0OffsetNs);
      Option(command, "--t0-policy", realT0Policy);
    } else {
      Option(command, "--fixed-t0-ns", realFixedT0Ns);
    }
  }

  std::cout << "\nCeeV0Analysis command:\n" << command.str() << "\n" << std::endl;
  if (dryRun) return 0;
  const int status = std::system(command.str().c_str());
  int result = 1;
  if (status == -1) std::cerr << "Could not launch CeeV0Analysis.\n";
  else if (WIFEXITED(status)) result = WEXITSTATUS(status);
  else if (WIFSIGNALED(status)) result = 128 + WTERMSIG(status);
  std::cout << "CeeV0Analysis exit code: " << result << std::endl;
  // The caller can use gSystem->Exit(runCeeV0Reco(...)) for batch-job exit status.
  return result;
}
