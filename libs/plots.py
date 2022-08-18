from typing import Dict
import matplotlib.pyplot as plt

def plot_plddt(plots_path: str, ranked_models_dict: Dict):
    
    plt.clf()
    for ranked, ranked_path in sorted(ranked_models_dict.items(), key=lambda x: int("".join([i for i in x[0] if i.isdigit()]))):
        plddt_list = []
        with open(ranked_path) as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        plt.plot(res_list, plddt_list, label=ranked)
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLLDT')
    plt.savefig(f'{plots_path}/pLDDT.png')

def create_plots(plots_path: str, ranked_models_dict: Dict):

    plot_plddt(plots_path=plots_path, ranked_models_dict=ranked_models_dict)