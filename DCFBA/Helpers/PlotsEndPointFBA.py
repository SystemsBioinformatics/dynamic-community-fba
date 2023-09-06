import matplotlib.pyplot as plt
from ..DynamicModels import EndPointFBA


def plot_biomasses(ep: EndPointFBA, ls: list[str] = []):
    FBAsol = ep.m_model.getSolutionVector(names=True)
    FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

    # Set model ids and initial biomasses
    modelIds = ep.m_model.get_model_ids()
    if len(ls) == 0:
        ls = modelIds

    for i, mid in enumerate(ls):
        ls = []
        ls.append(-1 * FBAsol[f"BM_{mid}_exchange"])
        for key, value in FBAsol.items():
            if key.startswith(f"BM_{mid}_time"):
                ls.append(value)
        ls.append(FBAsol[f"BM_{mid}_exchange_final"])
        plt.plot(
            list(range(0, len(ep.m_times) + 1)),
            ls,
            color=f"C{i}",
            label=f"{mid}",
        )

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()


def plot_metabolites(
    ep: EndPointFBA, speciesIdsConcentration: dict[str, float]
):
    FBAsol = ep.m_model.getSolutionVector(names=True)
    FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

    for i, sid in enumerate(speciesIdsConcentration.keys()):
        ls = []
        ls.append(speciesIdsConcentration[sid])
        for key, value in FBAsol.items():
            if key.startswith(f"{sid}_time"):
                ls.append(value)

        ls.append(FBAsol[f"{sid}_exchange_final"])
        plt.plot(
            list(range(0, len(ep.m_times) + 1)),
            ls,
            color=f"C{i}",
            label=f"{sid}",
        )

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()
