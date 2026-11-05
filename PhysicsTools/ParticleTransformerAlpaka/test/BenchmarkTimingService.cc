// Completion timestamps after event cleanup, for multi-stream wall throughput.
#include <chrono>
#include <fstream>
#include <mutex>
#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/ServiceRegistry/interface/StreamContext.h"
class BenchmarkTimingService {
  struct Completion { unsigned long long event; unsigned stream; long long ns; };
public:
  BenchmarkTimingService(edm::ParameterSet const& p, edm::ActivityRegistry& registry)
      : path_(p.getParameter<std::string>("output")) {
    registry.watchPostClearEvent(this,&BenchmarkTimingService::complete);
    registry.watchPostEndJob(this,&BenchmarkTimingService::write);
  }
  void complete(edm::StreamContext const& context) {
    auto now=std::chrono::steady_clock::now().time_since_epoch();
    auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();
    std::lock_guard<std::mutex> guard(mutex_);
    completions_.push_back({context.eventID().event(),context.streamID().value(),ns});
  }
  void write() {
    std::ofstream out(path_);
    for(auto const& x:completions_)
      out << "{\"event\":" << x.event << ",\"stream\":" << x.stream << ",\"ns\":" << x.ns << "}\n";
  }
private:
  std::string path_;
  std::mutex mutex_;
  std::vector<Completion> completions_;
};
DEFINE_FWK_SERVICE(BenchmarkTimingService);
