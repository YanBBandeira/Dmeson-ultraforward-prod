import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd 
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from pathlib import Path

plt.rcParams.update({    
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "figure.figsize": (5.7, 5.5),  # Ajustado ligeiramente para acomodar melhor os 2 painéis
    "figure.dpi": 300,
    "font.size": 15,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 12,
    "axes.linewidth": 2.0,
    
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.minor.width": 1.5,
    "ytick.minor.width": 1.5,
    "xtick.major.size": 7.0,
    "ytick.major.size": 7.0,
    "xtick.minor.size": 4.0,
    "ytick.minor.size": 4.0,
    
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.formatter.limits": [-2, 3],
    "axes.formatter.use_mathtext": True,
    "grid.linestyle": ":",
    "grid.linewidth": 0.6,
    "grid.alpha": 0.7,
    "grid.color": "gray",
    "lines.linewidth": 3.0,
    "legend.handlelength": 2.0,
    "legend.labelspacing": 0.05,
    "legend.borderaxespad": 0.1
})

BASE_DIR = Path(__file__).resolve().parents[2] if "__file__" in locals() else Path(".")

INPUT_PP_DIR = BASE_DIR / "Output" / "pp" / "Soft_pp" / "eCM_13000GeV"
INPUT_PP_FOTONPROD_DIR = BASE_DIR / "Output" / "pp" / "Soft_pp_fotonprod" / "eCM_13000GeV"
FIGURE_DIR = BASE_DIR /  "Analise" / "Figures"

FIGURE_DIR.mkdir(parents=True, exist_ok=True)


CHANNEL_STYLES = {
    "hadronico": {
        "color": "blue",
        "marker": "o",
        "linestyle": "-",
        "base_label": "Hadronic"
    },
    "foton": {
        "color": "red",
        "marker": "s",
        "linestyle": "-",
        "base_label": "Photon prod."
    }
}

RATIO_STYLE = {
    "color": "black",
    "marker": "D",
    "linestyle": "-",
    "label": "Ratio (Hadronic / Photon)"
}

datasets = {
    "D0": {
        "hadronico": INPUT_PP_DIR / "dNdeta_D0meson_pt.dat",
        "foton": INPUT_PP_FOTONPROD_DIR / "dNdeta_D0meson_pt.dat"
    },
    "Ds": {
        "hadronico": INPUT_PP_DIR / "dNdeta_Dsmeson_pt.dat",
        "foton": INPUT_PP_FOTONPROD_DIR / "dNdeta_Dsmeson_pt.dat"
    }
}

def carregardados(inputfile):
      if not inputfile.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {inputfile}")
      
      with open(inputfile) as f:
            linhas = [line.strip().split() for line in f
            if line.strip() and not line.startswith("#")]
            
      valx  = np.array([float(x[0]) for x in linhas])
      valdN = np.array([float(x[1]) for x in linhas])

      print(f"  -> Lido {inputfile.name}: {len(valx)} pontos | y: [{valx.min():.1f}, {valx.max():.1f}]")
                  
      return valx, valdN

def calcular_razao(x_h, dN_h, x_p, dN_p):
      """Função auxiliar para interpolar e calcular a razão Hadronic / Photon Prod."""
      idx_h = np.argsort(x_h)
      idx_p = np.argsort(x_p)
      x_h_sorted, dN_h_sorted = x_h[idx_h], dN_h[idx_h]
      x_p_sorted, dN_p_sorted = x_p[idx_p], dN_p[idx_p]

      mask = (x_h_sorted >= x_p_sorted.min()) & (x_h_sorted <= x_p_sorted.max())
      x_ratio = x_h_sorted[mask]
      dN_h_at = dN_h_sorted[mask]
      dN_p_at = np.interp(x_ratio, x_p_sorted, dN_p_sorted)

      valid = dN_p_at > 0
      x_ratio = x_ratio[valid]
      dN_h_at = dN_h_at[valid]
      dN_p_at = dN_p_at[valid]
      ratio = dN_h_at / dN_p_at
      return x_ratio, ratio

def gerar_plot_individual(meson_name, paths_dict, output_pdf_path):
      """Gera o gráfico individual (Hadrônico vs Fotonprodução) com painel de razão."""
      print(f"\nGerando gráfico individual para o méson {meson_name}...")
      
      x_h, dN_h = carregardados(paths_dict["hadronico"])
      x_p, dN_p = carregardados(paths_dict["foton"])
      
      fig = plt.figure()
      gs = GridSpec(2, 1, height_ratios=[3.6, 1], hspace=0)
      ax1 = fig.add_subplot(gs[0])
      ax2 = fig.add_subplot(gs[1])
      
      # Painel Superior
      for canal, filepath in paths_dict.items():
            st = CHANNEL_STYLES[canal]
            x, dN = carregardados(filepath)
            ax1.plot(x, dN, marker=st["marker"], color=st["color"], linestyle=st["linestyle"],
                     linewidth=2.0, markersize=5, markerfacecolor="none", markeredgewidth=1.2, 
                     alpha=0.85, label=f"Soft pp ({st['base_label']})", zorder=2)
            
      ax1.text(0.02, 0.97, rf"$\sqrt{{s}} = 13.0$ TeV, {meson_name}",
               horizontalalignment="left", verticalalignment="top", transform=ax1.transAxes)
      
      ax1.set_ylabel(r"$dN/dp_T$")
      ax1.set_yscale("log")
      ax1.set_xlim(5.0, 10.0)
      ax1.grid(linestyle=":", linewidth=0.6, alpha=0.7)
      ax1.set_xticklabels([])
      ax1.legend(frameon=False, fontsize=12, loc="best", labelspacing=0.1, handlelength=1.8, handletextpad=0.5)

      # Painel Inferior (Ratio)
      x_ratio, ratio = calcular_razao(x_h, dN_h, x_p, dN_p)
      ax2.plot(x_ratio, ratio, marker=RATIO_STYLE["marker"], color=RATIO_STYLE["color"], linestyle=RATIO_STYLE["linestyle"],
               linewidth=1.5, markersize=4, markerfacecolor="none", markeredgewidth=1.2, alpha=0.85, zorder=2)
      ax2.axhline(y=1.0, color="black", linestyle="--", linewidth=0.7, alpha=0.7)

      ax2.text(0.02, 0.92, "Hadronic / Photon Prod.", verticalalignment="top", horizontalalignment="left", transform=ax2.transAxes, fontsize=11)
      ax2.set_xlabel(r"$p_T [\mathrm{GeV}]$")
      ax2.set_ylabel(r"Ratio")
      ax2.set_ylim(0, 2)      
      ax2.set_xlim(5.0, 10.0)
      ax2.grid(linestyle=":", linewidth=0.6, alpha=0.7)

      plt.tight_layout()
      with PdfPages(output_pdf_path) as pp:
            pp.savefig(fig, bbox_inches="tight", dpi=300)
      plt.close(fig)
      print(f" -> Salvo em: {output_pdf_path}")

def gerar_plot_conjunto_com_razao(datasets, output_pdf_path):
      """Gera o gráfico conjunto reunindo D0 e Ds com painel superior (distribuições) e inferior (razões)."""
      print(f"\nGerando gráfico conjunto (D0 e Ds) com painel de razão...")
      
      fig = plt.figure()
      gs = GridSpec(2, 1, height_ratios=[3.6, 1], hspace=0)
      ax1 = fig.add_subplot(gs[0])
      ax2 = fig.add_subplot(gs[1])
      
      # Variações de marcadores/cores customizadas para o gráfico conjunto para distinguir D0 de Ds mantendo o padrão
      conjunto_styles = {
          ("D0", "hadronico"): {"color": "blue", "marker": "o", "label": r"$D^0$ (Hadronic)"},
          ("D0", "foton"):     {"color": "cyan", "marker": "^", "label": r"$D^0$ (Photon prod.)"},
          ("Ds", "hadronico"): {"color": "red", "marker": "s", "label": r"$D_s$ (Hadronic)"},
          ("Ds", "foton"):     {"color": "magenta", "marker": "v", "label": r"$D_s$ (Photon prod.)"}
      }
      
      ratio_styles_conjunto = {
          "D0": {"color": "blue", "marker": "D", "label": r"$D^0$ Ratio"},
          "Ds": {"color": "red", "marker": "s", "label": r"$D_s$ Ratio"}
      }
      
      # 1. Plotagem do Painel Superior (Distribuições)
      for meson, paths in datasets.items():
          for canal, filepath in paths.items():
              if filepath.exists():
                  x, dN = carregardados(filepath)
                  cfg = conjunto_styles[(meson, canal)]
                  ax1.plot(x, dN, marker=cfg["marker"], color=cfg["color"], linestyle="-",
                           linewidth=2.0, markersize=5, markerfacecolor="none", markeredgewidth=1.2,
                           alpha=0.85, label=cfg["label"])
      ax1.text(0.02, 0.97, r"$\sqrt{s} = 13.0$ $\, \mathrm{TeV}$, $D^0$ \& $D_s$",
              horizontalalignment="left", verticalalignment="top", transform=ax1.transAxes)
      
      ax1.set_ylabel(r"$dN/dpt$")
      ax1.set_yscale("log")
      ax1.set_xlim(0.0, 40.0)
      ax1.grid(linestyle=":", linewidth=0.6, alpha=0.7)
      ax1.set_xticklabels([])
      ax1.legend(frameon=False, fontsize=10, loc="best", labelspacing=0.1, handlelength=1.5, ncol=2)

      # 2. Plotagem do Painel Inferior (Razões Hadronico / Fotonprod para D0 e Ds)
      for meson, paths in datasets.items():
          if paths["hadronico"].exists() and paths["foton"].exists():
              x_h, dN_h = carregardados(paths["hadronico"])
              x_p, dN_p = carregardados(paths["foton"])
              
              x_ratio, ratio = calcular_razao(x_h, dN_h, x_p, dN_p)
              rcfg = ratio_styles_conjunto[meson]
              
              ax2.plot(x_ratio, ratio, marker=rcfg["marker"], color=rcfg["color"], linestyle="-",
                       linewidth=1.5, markersize=4, markerfacecolor="none", markeredgewidth=1.2,
                       alpha=0.85, label=rcfg["label"])

      ax2.axhline(y=1.0, color="black", linestyle="--", linewidth=0.7, alpha=0.7)
      ax2.set_xlabel(r"$p_T [\mathrm{GeV}]$")
      ax2.set_ylabel(r"Ratio")
      ax2.set_ylim(0, 2)      
      ax2.set_xlim(5.0, 10.0)
      ax2.grid(linestyle=":", linewidth=0.6, alpha=0.7)
      ax2.legend(frameon=False, fontsize=10, loc="best", labelspacing=0.1, handlelength=1.5)

      plt.tight_layout()
      with PdfPages(output_pdf_path) as pp:
            pp.savefig(fig, bbox_inches="tight", dpi=300)
      plt.close(fig)
      print(f" -> Salvo em: {output_pdf_path}")

if __name__ == "__main__":
    # 1. Gráfico individual para o D0 (com razão)
    gerar_plot_individual("D0", datasets["D0"], FIGURE_DIR / "plot_dNpt_D0_comparison.pdf")
    
    # 2. Gráfico individual para o Ds (com razão)
    gerar_plot_individual("Ds", datasets["Ds"], FIGURE_DIR / "plot_dNdpt_Ds_comparison.pdf")
    
    # 3. Gráfico conjunto (D0 e Ds juntos) com painel de razão integrado
    gerar_plot_conjunto_com_razao(datasets, FIGURE_DIR / "plot_dNdpt_D0_Ds_combined_ratio.pdf")
    
    print("\nTodos os gráficos foram gerados e salvos com sucesso na pasta Figures!")