import json
from collections import OrderedDict
from matplotlib import pyplot as plt

njets = 5446731

use_jets = False  # if false, dividing by total orbits
logy = True
# modules = ["SoftJetIdSC4"]
modules = {
    "scPhase2PFRawToDigiAlpaka": "Puppi unpacking",
    "scPhase2VertexRawToDigiAlpaka": "Vertex unpacking",
    "scPhase2SC4PFAlpaka": "SC4 jet clustering",
    "SoftJetIdSC4": "Inference",
    "cleanup": "Cleanup",
    "other": "Other",
    # "total": "Total time"
}

colors = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
colors = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
            "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
styles = ["-", "--", "-.", ":"]

labels_legends = {
    # "part": "ParT",
    # "test": "ParT INT8",
    #"graphpart": "Graph ParT",
    # "deepsets_optim": "Deepsets",
    # "graphdeepsets": "Graph Deepsets",
    #"partJJ48": "Part 48 const",
    #"partJJ12": "Part 12 const, 16 components",
    #"part_2p5k": "Part 2.5k",
    # "part_2p5k_throughput": "Part 2.5k, throughput opt.",
    #"partJJ12_48const": "Part 12 const, 48 components",

    # new models
    "deepsets_2p5k": "Deepsets 2.5k (NEW)",
    "deepsets_10k": "Deepsets 10k",
    "deepsets_50k": "Deepsets 50k",
    "deepsets_300k": "Deepsets 300k",
    "deepsets_2p2m": "Deepsets 2.2m",
    "pnet_2p5k": "PNet 2.5k",
    "pnet_10k": "PNet 10k",
    "pnet_50k": "PNet 50k",
    "pnet_300k": "PNet 300k",
    "pnet_2p2m": "PNet 2.2m",
    "part_2p5k": "Part 2.5k",
    "part_50k": "Part 50k",
    "part_10k": "Part 10k",
    "part_300k": "Part 300k",
    "part_2p2m": "Part 2.2m",
}

def get_total_time(modules, d):
    if type(modules) == list or type(modules) == dict:
        return sum([get_total_time(elem, d) for elem in modules])
    else:
        if not modules == "total":
            inference_index = [elem["label"] for elem in d["modules"]].index(modules)
            time = d["modules"][inference_index]["time_real"]
            nevents = d["modules"][inference_index]["events"]
        else:
            time = d["total"]["time_real"]
            nevents = d["total"]["events"]
        if not use_jets:
            return time / nevents
        else:
            return time / njets


def create_dict(keys=["ParT", "Deepsets"]):
    return {
        key: OrderedDict()
        for key in keys
    }

def plot(labels_results, x_title, plotname, inner_text, legends=None):
    fig, ax = plt.subplots()#layout='constrained')
    for i, (label, results) in enumerate(labels_results.items()):
        if legends:
            label = legends[i]
        color = colors[i % len(colors)]
        style = styles[i // len(colors)]
        plt.plot(results.keys(), results.values(), label=label, color=color, linestyle=style)
    ax.set_xlabel(x_title)
    if not use_jets:
        suf = ""
        ax.set_ylabel("Time per orbit [ms]")
    else:
        suf = "_jets"
        ax.set_ylabel("Time per jet [ms]")
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.text(0.6, 0.85, inner_text, transform=ax.transAxes)
    plt.savefig(f"plots/{plotname}{suf}.pdf")

    fig, ax = plt.subplots()#layout='constrained')
    for i, (label, results) in enumerate(labels_results.items()):
        if legends:
            label = legends[i]
    
        color = colors[i % len(colors)]
        style = styles[i // len(colors)]
        plt.plot(results.keys(), [1000 / elem for elem in results.values()], label=label, color=color, linestyle=style)
    ax.set_xlabel(x_title)
    if not use_jets:
        ax.set_ylabel("Orbits / s")
    else:
        ax.set_ylabel("Jets / s")
    if logy:
        ax.set_yscale("log")
    ax.text(0.6, 0.85, inner_text, transform=ax.transAxes)
    ax.legend(fontsize=8)
    plt.savefig(f"plots/{plotname}{suf}_freq.pdf")
    plt.close(fig)
    # fig, ax = plt.subplots()#layout='constrained')
    # plt.plot(results_thread.keys(), results_thread.values())
    # ax.set_xlabel(x_title)
    # ax.set_ylabel("Time_thread per orbit [ms]")
    # ax.text(0.6, 0.85, inner_text, transform=ax.transAxes)
    # plt.savefig(f"{plotname}_thread.pdf")



def plot_modules(labels_results, x_title, plotname, inner_text, legends=None):

    fig, ax = plt.subplots()
    total_bins = len(labels_results)
    total_modules = len(next(iter(labels_results.values())))
    total_points = total_bins * (1 + total_modules)

    for i, (label, results) in enumerate(labels_results.items()):
        for imodule, (module_name, result) in enumerate(results.items()):
            if legends:
                module_name = legends[imodule]
            if i == 0:
                plt.bar(
                    x=i * (total_modules + 1) + imodule,
                    height=result,
                    width=0.9,
                    color=colors[imodule],
                    label=module_name
                )
            else:
                plt.bar(
                    x=i * (total_modules + 1) + imodule,
                    height=result,
                    width=0.9,
                    color=colors[imodule],
                )
    
    ## Set Axes Labels
    plt.xlabel(x_title)
    plt.ylabel("Time per orbit [ms]")

    plt.xticks(
        ticks=[total_modules / 2. - 0.5 + (total_modules + 1) * i for i in range(total_bins)],
        labels=labels_results.keys()
    )
    plt.legend(loc='best')
    ax.text(0.6, 0.85, inner_text, transform=ax.transAxes)
    if logy:
        ax.set_yscale("log")
    plt.savefig(f"plots/{plotname}.pdf")
    plt.close(fig)

results_labels = create_dict(labels_legends.keys())

for label in results_labels:
    for backend in ["cuda_async"]: # serial_sync #
        for numThreads in [1]: #1 #2 4 8 16
            for numFwkStreams in [1]: #1 #2 4 8 16
                for batchSize in [32, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,]:# 65536, 262144, 1000000, 2000000, 4000000, 8000000, 16000000]:
                # for batchSize in [32, 512, 8192, 32768, 262144, 524288, 1048576, 2097152]:
                    try:
                        with open(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json") as f:
                            d = json.load(f)
                        key = f"{batchSize}"
                        results_labels[label][key] = get_total_time(modules, d)

                    except FileNotFoundError:
                        print(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json not found")
                    # results_thread[key] = d["modules"][6]["time_thread"] / d["modules"][6]["events"]

plot(results_labels, "Batch size", "results_cuda_batch", "cuda_async \n1 thread \n1 stream", list(labels_legends.values()))

results_labels = create_dict(labels_legends.keys())

batch_sizes = [1024]
for label in results_labels:
    for backend in ["cuda_async"]: # serial_sync #
        for numThreads in [1, 2, 4, 8]: #1 #2 4 8 16
            for numFwkStreams in [1, 2, 4, 8]: #1 #2 4 8 16
                for batchSize in batch_sizes:
                    try:
                        with open(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json") as f:
                            d = json.load(f)
                        key = f"({numThreads}, {numFwkStreams})"
                        results_labels[label][key] = get_total_time(modules, d)
                    except FileNotFoundError:
                        print(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json not found")

label = "cuda_async"
if len(batch_sizes) == 0:
    label += f"\nbatch size = {batch_sizes[0]}"
plot(results_labels, "(# threads, # streams)", "results_cuda_threads", "cuda_async \nbatch size = 32", list(labels_legends.values()))


results_labels = create_dict(labels_legends.keys())

for label in results_labels:
    for backend in ["serial_sync", "cuda_async"]: # serial_sync #
        for numThreads in [1, 8]: #1 #2 4 8 16
            for numFwkStreams in [1, 8]: #1 #2 4 8 16
                for batchSize in [32]:
                    try:
                        if "cuda" in backend and ((numThreads != 1) or (numFwkStreams != 1)):
                            continue
                        with open(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json") as f:
                            d = json.load(f)
                        key = f"({backend}, {numThreads}, {numFwkStreams})"
                        results_labels[label][key] = get_total_time(modules, d)
                        # results_thread[key] = d["modules"][6]["time_thread"] / d["modules"][6]["events"]

                    except FileNotFoundError:
                        print(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json not found")


plot(results_labels, "(backend, # threads, # streams)", "results_cpu_gpu", "batch size = 32", list(labels_legends.values()))


results_labels = create_dict(labels_legends.keys())

for label in results_labels:
    for backend in ["cuda_async"]: # serial_sync #
        for numThreads in [1]: #1 #2 4 8 16
            for numFwkStreams in [1]: #1 #2 4 8 16
                # for batchSize in [32, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]:
                # for batchSize in [32, 512, 8192, 32768, 262144, 524288, 1048576, 2097152]:
                for batchSize in [32, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,]:
                    try:
                        if "cuda" in backend and ((numThreads != 1) or (numFwkStreams != 1)):
                            continue
                        with open(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json") as f:
                            d = json.load(f)
                        key = f"{batchSize}"
                        results_labels[label][key] = get_total_time(modules, d)
                        # results_thread[key] = d["modules"][6]["time_thread"] / d["modules"][6]["events"]

                    except FileNotFoundError:
                        print(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json not found")

plot(results_labels, "Batch size", "results_cmssw", "cuda_async \n1 thread \n1 stream", list(labels_legends.values()))


for label, legend in labels_legends.items():
    results = {}
    for backend in ["cuda_async"]: # serial_sync #
        for numThreads in [1]: #1 #2 4 8 16
            for numFwkStreams in [1]: #1 #2 4 8 16
                # for batchSize in [32, 256, 1024, 4096, 16384, 65536, 262144]:
                for batchSize in [32, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,]:
                # for batchSize in [32, 512, 8192, 32768, 262144, 524288, 1048576, 2097152]:
                    try:
                        if "cuda" in backend and ((numThreads != 1) or (numFwkStreams != 1)):
                            continue
                        with open(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json") as f:
                            d = json.load(f)
                        key = f"{batchSize}"
                        results[key] = {}
                        for module in modules:
                            results[key][module] = get_total_time(module, d)
                        # results_thread[key] = d["modules"][6]["time_thread"] / d["modules"][6]["events"]

                    except FileNotFoundError:
                        print(f"results_{label}/resources.cms-p2scout-ngt-01.{backend}.{numThreads}.{numFwkStreams}.{batchSize}.json not found")
    path_label = "_".join(label.split("/"))
    plot_modules(results, "Batch size", f"results_modules_{path_label}", f"cuda_async \n1 thread, 1 stream \n{legend}", list(modules.values()))
