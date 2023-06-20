from cbmpy.CBModel import Model, Reaction, Species, Reagent
from KineticModel import KineticModel
import cbmpy
import sys
import pandas as pd
import numpy as np
import math

# sys.path is a list of absolute path strings
sys.path.append("code/tests/helpers")  # <-- relative path
import model_a
import model_b
import model_c

import matplotlib.pyplot as plt


def dynamic_fba(
    kinetic_model: KineticModel,
    biomass_reaction_id: str,
    time_steps: list[float],
    initial_conditions: dict[str, float],
):
    if "biomass" not in initial_conditions.keys():
        raise Exception("Biomass must be defined in the initial condition")

    model = kinetic_model.get_model()
    update_conditions = initial_conditions

    # If we have Vmax of a reaction set this to be the upper bound
    # Else just use the model defined upper bound
    for rid, tt in kinetic_model.get_model_kinetics().items():
        reaction: Reaction = model.getReaction(rid)
        reaction.setUpperBound(tt[0])

    species_reaction_dict: dict[
        Species, list[Reaction]
    ] = create_species_reactions_dict(model)

    for species in model.species:
        if species.id not in update_conditions.keys():
            update_conditions[species.id] = 0

    t_old = 0

    for t in time_steps:
        dt = t - t_old
        Xt = update_conditions["biomass"]
        reaction: Reaction
        for reaction in model.reactions:
            if reaction.id not in model.getExchangeReactionIds():
                reagents: list[Reagent] = reaction.reagents
                slowest_reaction = 1e10
                for reagent in reagents:
                    species_id: str = reagent.getSpecies()
                    species = model.getSpecies(species_id)
                    if species.getCompartmentId() == "e":
                        reaction_upper_bound = calc_reaction_upper_bound(
                            kinetic_model,
                            reaction.id,
                            reaction.getUpperBound(),
                            update_conditions[species.id],
                            Xt,
                            dt,
                        )
                        if reaction_upper_bound < slowest_reaction:
                            reaction.setUpperBound(reaction_upper_bound)
                            slowest_reaction = reaction_upper_bound
        print("-----------")
        print(model.getReaction("R_GLCpts").getUpperBound())
        print("-----------")

        # TODO speed this up
        # for species_id, ls in species_reaction_dict.items():
        #     for reaction in ls:
        #         reaction_upper_bound = calc_reaction_upper_bound(
        #             kinetic_model,
        #             reaction.id,
        #             reaction.getUpperBound(),
        #             update_conditions[species_id],
        #             Xt,
        #             dt,
        #         )
        #         reaction.setUpperBound(reaction_upper_bound)

        solution = cbmpy.doFBA(model)
        if math.isnan(solution):
            print(t)
            break

        # TODO fix this to dict maybe?
        FBAsol = model.getSolutionVector(names=True)
        FBAsol = pd.DataFrame(
            zip(FBAsol[1], FBAsol[0]), columns=["ID", "flux"]
        )
        print(FBAsol)
        # calculate new concentrations
        # TODO make function out of this
        for species_id, ls in species_reaction_dict.items():
            new_concentration = update_conditions[species_id]
            reaction: Reaction
            for reaction_id in ls:
                new_concentration = new_concentration - (
                    float(FBAsol.loc[FBAsol["ID"] == reaction_id, "flux"])
                    * Xt
                    * dt
                )
            # update in dictonary
            update_conditions[species_id] = new_concentration
        # Update Biomass
        update_conditions["biomass"] = Xt + (
            float(FBAsol.loc[FBAsol["ID"] == biomass_reaction_id]["flux"])
            * Xt
            * dt
        )
        t_old = t


def create_species_reactions_dict(
    model: Model,
) -> dict[Species, list[Reaction]]:
    """We need to access the species and the reactions that update the species
    a lot, so we start by building the species reactions dict
    which we can than easily access in the dynamic fba function

    Args:
        model (Model): _description_
    """
    dict = {}
    species: Species
    for species in model.species:
        if species.getCompartmentId() == "e":
            dict[species.id] = []
            for reaction in species.isReagentOf():
                if reaction not in model.getExchangeReactionIds():
                    dict[species.id].append(reaction)
    return dict


def calc_reaction_upper_bound(kinetic_model: KineticModel, rid, ub, S, X, dt):
    km = kinetic_model.get_reaction_km(rid)
    if km is None:
        # Calculate Without kinetics
        v_hat = S / (X * dt)
        if ub <= v_hat:
            return ub
        return v_hat
    else:
        print("2222222")
        print(ub, S, km, S)
        print(ub * (S / (km + S)))
        print("2222222")

        return ub * (S / (km + S))


def dynamic_fba_fault():
    model: Model = model_c.build_model_C()
    # Imprt S into the cell this is positive
    import_S: Reaction = model.getReaction("R_1")

    # Lets see if we can set concentrations by setting the exchange
    s_exchange: Reaction = model.getReaction("S_exchange")

    # Inirial concentration 5?

    bm_exchange: Reaction = model.getReaction("BM_e_C_exchange")

    bm_exchange.setLowerBound(0)
    # No exchanges?
    bm_exchange.setUpperBound(1000)
    concentration_S_t = 10
    concentration_BM_t = 0.1
    s_exchange.setLowerBound(-concentration_S_t)
    # No exchanges?
    s_exchange.setUpperBound(5)
    y1 = []
    y2 = []
    ts = []
    prev_t = 0

    for t in np.linspace(0, 15, 500):
        dt = t - prev_t
        v_s = 10 * (concentration_S_t / (5 + concentration_S_t))
        import_S.setUpperBound(v_s)

        solution = cbmpy.doFBA(model)
        if math.isnan(solution):
            print(t)
            print(concentration_S_t)
            print(concentration_BM_t)
            break
        FBAsol = model.getSolutionVector(names=True)
        FBAsol = pd.DataFrame(
            zip(FBAsol[1], FBAsol[0]), columns=["ID", "flux"]
        )

        # calculate new concentrations
        concentration_S_t = concentration_S_t - (
            float(FBAsol.loc[FBAsol["ID"] == "R_1", "flux"])
            * concentration_BM_t
            * dt
        )
        concentration_BM_t = concentration_BM_t + (
            float(FBAsol.loc[FBAsol["ID"] == "BM_e_C_exchange"]["flux"])
            * concentration_BM_t
            * dt
        )
        y1.append(concentration_S_t)
        y2.append(concentration_BM_t)
        ts.append(t)

        prev_t = t
    return [ts, y1, y2]


# ts, y1, y2 = dynamic_fba()
# print(len(ts))
# print(len(y1))
# print(len(y2))
# ax = plt.subplot(111)
# ax.plot(ts, y2)
# ax2 = plt.twinx(ax)
# ax2.plot(ts, y1, color="r")

# ax.set_ylabel("Biomass", color="b")
# ax2.set_ylabel("Glucose", color="r")

# plt.show()


def dynamic_ecoli_core():
    model: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")

    glucose: Reaction = model.getReaction("R_GLCpts")
    concentration_glucose_t = 10
    concentration_BM_t = 0.1
    # bm = model.getReaction("R_BIOMASS_Ecoli_core_w_GAM")
    # aa: Objective = model.getActiveObjective()

    glucose_max_import = (
        10 * concentration_glucose_t / (5 + concentration_glucose_t)
    )
    glucose.setUpperBound(glucose_max_import)

    y1 = []
    y2 = []
    ts = []
    t_old = 0
    bm_fluxes = []
    for t in np.linspace(0, 15, 100):
        dt = t - t_old
        print(dt)
        v_s = 10 * (concentration_glucose_t / (5 + concentration_glucose_t))
        glucose.setUpperBound(v_s)
        # v_hat = (
        #     glucose.getLowerBound()
        #     * (concentration_glucose_t)
        #     / (5 + concentration_glucose_t)
        # )
        # if v_hat <= v_s:
        #     print(concentration_glucose_t)
        #     print("----")
        #     print(t)
        #     glucose.setLowerBound(v_hat)
        solution = cbmpy.doFBA(model)
        if (
            math.isnan(solution)
            or solution <= 0.00001
            or concentration_glucose_t < 0.00001
        ):
            print(t)
            print(concentration_glucose_t)
            print(concentration_BM_t)
            break
        FBAsol = model.getSolutionVector(names=True)
        FBAsol = pd.DataFrame(
            zip(FBAsol[1], FBAsol[0]), columns=["ID", "flux"]
        )
        print(FBAsol)
        # calculate new concentrations
        concentration_glucose_t = concentration_glucose_t - (
            float(FBAsol.loc[FBAsol["ID"] == "R_GLCpts", "flux"])
            * concentration_BM_t
            * dt
        )
        print(concentration_glucose_t)
        if t > 0.3:
            raise Exception
        concentration_BM_t = concentration_BM_t + (
            float(
                FBAsol.loc[FBAsol["ID"] == "R_BIOMASS_Ecoli_core_w_GAM"][
                    "flux"
                ]
            )
            * concentration_BM_t
            * dt
        )
        bm_fluxes.append(solution)
        y1.append(concentration_glucose_t)
        y2.append(concentration_BM_t)
        ts.append(t)
        t_old = t
    return [ts, y1, y2, bm_fluxes]


# ts, y1, y2, bm_fluxes = dynamic_ecoli_core()

# print(ts)
# print((y1))
# print((y2))
# ax = plt.subplot(111)
# ax.plot(ts, y2)
# ax2 = plt.twinx(ax)
# ax2.plot(ts, y1, color="r")

# ax.set_ylabel("Biomass", color="b")
# ax2.set_ylabel("Glucose", color="r")


# plt.show()


model: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
ex_r: Reaction = model.getReaction("R_EX_glc__D_e")
ex_r.is_exchange = True
kinetic_model: KineticModel = KineticModel(model, {"R_GLCpts": [10, 5]})

dynamic_fba(
    kinetic_model,
    "R_BIOMASS_Ecoli_core_w_GAM",
    np.linspace(0, 15, 100),
    {"M_glc__D_e": 10, "biomass": 0.1},
)
