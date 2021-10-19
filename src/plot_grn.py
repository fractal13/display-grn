#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import matplotlib as mplt

def draw_interaction(ax, ar):
    center = (0.5, 0.5)
    a_color = (43/255,102/255,158/255)
    r_color = (235/255,117/255,110/255)
    t_color = (100/255,100/255,100/255)
    
    a = ar[0]
    r = ar[1]

    scale = 0.45

    if a > r:
        circle = mplt.patches.Circle(center, a*scale, color=a_color)
        ax.add_patch(circle)
        circle = mplt.patches.Circle(center, r*scale, color=r_color)
        ax.add_patch(circle)
    elif r > a:
        circle = mplt.patches.Circle(center, r*scale, color=r_color)
        ax.add_patch(circle)
        circle = mplt.patches.Circle(center, a*scale, color=a_color)
        ax.add_patch(circle)
    else:
        if r >= 1.0:
            pass
        else:
            circle = mplt.patches.Circle(center, r*scale, color=t_color)
            ax.add_patch(circle)
        
    return

g_expected_genes = ["cad", "bcd", "tor", "nos", "gt", "tll", "tsl", "hb", "kr", "kni", "ftz", "eve", "run", "h", "prd"]
g_gene_map = {
    "cad": "cad",
    "Cad": "cad",
    "bcd": "bcd",
    "Bcd": "bcd",
    "tor": "tor",
    "nos": "nos",
    "Nos": "nos",
    "gt": "gt",
    "Gt": "gt",
    "tll": "tll",
    "Tsl": "tsl",
    "hb": "hb",
    "Hb": "hb",
    "Kr": "kr",
    "kr": "kr",
    "kni": "kni",
    "Kni": "kni",
    "ftz": "ftz",
    "ftz": "ftz",
    "eve": "eve",
    "Eve": "eve",
    "run": "run",
    "Run": "run",
    "h": "h",
    "Prd": "prd",
}

def update_genes_and_network(genes, network):
    genes = [ g_gene_map[gene] for gene in genes ]
    locations = {}
    for gene in g_expected_genes:
        if gene in genes:
            i = genes.index(gene)
        else:
            i = -1
        locations[gene] = i
    for gene in genes:
        if gene not in g_expected_genes:
            raise Exception("gene ({}) not expected".format(gene))


    updated_network = []
    for gene1 in g_expected_genes:
        row = []
        i1 = locations[gene1]
        for gene2 in g_expected_genes:
            i2 = locations[gene2]
            if i1 < 0 or i2 < 0:
                a, r = 1.0, 1.0
            else:
                a, r = network[i1][i2]
            row.append((a,r))
        updated_network.append(row)
    updated_genes = g_expected_genes[:]
    return updated_genes, updated_network

def plot_grn(genes, network, file_name):
    genes, network = update_genes_and_network(genes, network)

    fig_width = 6
    fig_height = 6
    n_genes = len(genes)

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='white')
    ax_tmp = fig.subplots(nrows=n_genes, ncols=n_genes, sharex='all', sharey='all', subplot_kw={'frameon':True})

    all_ax = []
    for row in range(n_genes):
        for col in range(n_genes):
            ax_tmp[row][col].xaxis.set_visible(False)
            ax_tmp[row][col].yaxis.set_ticks([])
            all_ax.append(ax_tmp[row][col])
            if row == 0:
                all_ax[len(all_ax)-1].set_title('{}'.format(genes[col]))
            if col == 0:
                all_ax[len(all_ax)-1].set_ylabel('{}'.format(genes[row]))

    fig.suptitle("Regulated")
    fig.supylabel("Regulator")
    position = 0
    for row in range(n_genes):
        for col in range(n_genes):
            ar = network[row][col]
            ax = all_ax[position]
            draw_interaction(ax, ar)
            position += 1

    fig.savefig(file_name,  dpi=100)
    return


def main2():
    genes = ["a", "b", "c"]
    network = [[(1,0), (1,1), (0,1)],
               [(0.5,0.5), (0.8,0.3), (0.8,0.95)],
               [(0,0.8), (0.2,1), (1,0.2)],
    ]
    plot_grn(genes, network, "foo.png")
    return

def Jaeger():
    genes = ["hb", "bcd", "gt", "kr", "tll", "kni"]
    network = [
        # hb        bcd        gt         kr         tll        kni
        [(1.0,0.0), (1.0,1.0), (0.0,0.9), (0.7,0.8), (1.0,1.0), (0.0,1.0) ], # hb
        [(1.0,0.0), (1.0,1.0), (1.0,0.0), (1.0,0.0), (1.0,1.0), (1.0,0.0) ], # bcd
        [(1.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (1.0,1.0), (0.0,0.7) ], # gt
        [(0.0,0.7), (1.0,1.0), (0.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,1.0) ], # kr
        [(1.0,1.0), (1.0,1.0), (0.0,0.7), (1.0,1.0), (1.0,1.0), (0.0,0.7) ], # tll
        [(0.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,0.7), (1.0,1.0), (1.0,1.0) ], # kni
    ]
    plot_grn(genes, network, "Jaeger.png")
    return

def Tanevski():
    genes = ["cad", "bcd", "nos", "gt", "tll", "hb", "kr", "kni", "ftz", "eve", "run", "h"]
    network = [
        # cad       bcd        nos        gt         tll        hb         kr         kni        ftz        eve        run        h
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # cad
        [(0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,0.0), (0.8,0.9), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # bcd
        [(1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # nos
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # gt
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0) ], # tll
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,0.0), (1.0,1.0) ], # hb
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0) ], # kr
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0) ], # kni
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # ftz
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (0.0,1.0), (1.0,0.0) ], # eve
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (0.0,1.0) ], # run
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0) ], # h
    ]
    plot_grn(genes, network, "Tanevski.png")
    return

def Rodriguez():
    genes = ["cad", "bcd", "nos", "gt", "tll", "hb", "kr", "kni"]
    network = [
        # cad       bcd        nos        gt         tll        hb         kr         kni        
        [(1.0,0.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (0.0,1.0), (1.0,0.0), (1.0,0.0) ], # cad
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,0.0), (1.0,0.0), (1.0,0.0) ], # bcd
        [(1.0,0.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0) ], # nos
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (0.0,1.0) ], # gt
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (0.0,1.0) ], # tll
        [(0.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (0.0,1.0) ], # hb
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,0.0), (0.0,1.0) ], # kr
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (0.0,1.0), (1.0,0.0) ], # kni
    ]
    plot_grn(genes, network, "Rodriguez.png")
    return

def Dilao():
    genes = ["cad", "bcd", "tor", "gt", "tll", "hb", "kr", "kni"]
    network = [
        # cad       bcd        tor        gt         tll        hb         kr         kni        
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0) ], # cad
        [(0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (1.0,0.0), (1.0,0.0), (1.0,0.0) ], # bcd
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (1.0,1.0), (1.0,1.0) ], # tor
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (0.0,1.0) ], # gt
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (0.0,1.0), (0.0,1.0) ], # tll
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (0.0,1.0) ], # hb
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0) ], # kr
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0) ], # kni
    ]
    plot_grn(genes, network, "Diloa.png")
    return

def Levine():
    genes = ["cad", "bcd", "tor", "gt", "tll", "hb", "kr", "kni"]
    network = [
        # cad       bcd        tor        gt         tll        hb         kr         kni        
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # cad
        [(0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,0.0), (1.0,1.0), (1.0,0.0) ], # bcd
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # tor
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0) ], # gt
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0) ], # tll
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,0.0), (0.8,0.6), (0.0,1.0) ], # hb
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0) ], # kr
        [(1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (1.0,1.0), (0.0,1.0), (1.0,1.0), (1.0,1.0) ], # kni
    ]
    plot_grn(genes, network, "Levine.png")
    return

def main():
    Jaeger()
    Tanevski()
    Rodriguez()
    Dilao()
    Levine()
    return

if __name__ == "__main__":
    main()
