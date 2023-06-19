from cbmpy.CBModel import Model, Reaction, Reagent
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


def dynamic_fba():
    model: Model = model_c.build_model_C()
    # Imprt S into the cell this is positive
    import_S: Reaction = model.getReaction("R_1")

    # Lets see if we can set concentrations by setting the exchange
    s_exchange: Reaction = model.getReaction("S_exchange")

    # Inirial concentration 5?
    s_exchange.setLowerBound(-10)
    # No exchanges?
    s_exchange.setUpperBound(5)

    bm_exchange: Reaction = model.getReaction("BM_e_C_exchange")

    bm_exchange.setLowerBound(0)
    # No exchanges?
    bm_exchange.setUpperBound(1000)
    concentration_S_t = 10
    concentration_BM_t = 0.1

    y1 = []
    y2 = []
    ts = []
    prev_t = 0

    # For every hour?
    # [x * 0.01 for x in range(1, 10000)]:
    for t in np.linspace(0, 15, 100):
        dt = t - prev_t
        v_s = 10 * (concentration_S_t / (5 + concentration_S_t))
        print(v_s)
        import_S.setUpperBound(v_s)
        # v_hat = calculate_v_hat(concentration_S_t, concentration_BM_t, dt)
        # if v_hat <= import_S.getUpperBound():
        #     print("----")
        #     print(t)
        #     import_S.setUpperBound(v_hat)
        solution = cbmpy.doFBA(model)
        if solution == 0:
            print(t)
            print(concentration_S_t)
            print(concentration_BM_t)
            break
        FBAsol = model.getSolutionVector(names=True)
        FBAsol = pd.DataFrame(
            zip(FBAsol[1], FBAsol[0]), columns=["ID", "flux"]
        )
        print(FBAsol)
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


def calculate_v_hat(St, Xt, dt):
    return St / (Xt * dt)


ts, y1, y2 = dynamic_fba()
print(len(ts))
print(len(y1))
print(len(y2))
ax = plt.subplot(111)
ax.plot(ts, y2)
ax2 = plt.twinx(ax)
ax2.plot(ts, y1, color="r")

ax.set_ylabel("Biomass", color="b")
ax2.set_ylabel("Glucose", color="r")

plt.show()


def dynamic_ecoli_core():
    model: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")

    glucose: Reaction = model.getReaction("R_EX_glc__D_e")
    concentration_glucose_t = 10
    concentration_BM_t = 0.1

    glucose_max_import = (
        -10 * concentration_glucose_t / (5 + concentration_glucose_t)
    )
    glucose.setLowerBound(glucose_max_import)

    y1 = []
    y2 = []
    ts = []
    t_min_1 = 0
    for t in np.linspace(0, 15, 100):
        print(t)
        v_s = -10 * (concentration_glucose_t / (5 + concentration_glucose_t))
        glucose.setLowerBound(v_s)
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
        if solution <= 0.00001 or concentration_glucose_t < 0.00001:
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
        concentration_glucose_t = concentration_glucose_t + (
            float(FBAsol.loc[FBAsol["ID"] == "R_EX_glc__D_e", "flux"])
            * concentration_BM_t
            * (t - t_min_1)
        )
        concentration_BM_t = concentration_BM_t + (
            float(
                FBAsol.loc[FBAsol["ID"] == "R_BIOMASS_Ecoli_core_w_GAM"][
                    "flux"
                ]
            )
            * concentration_BM_t
            * (t - t_min_1)
        )

        y1.append(concentration_glucose_t)
        y2.append(concentration_BM_t)
        ts.append(t)
        t_min_1 = t

    return [ts, y1, y2]


# ts, y1, y2 = dynamic_ecoli_core()

# # print(ts)
# print((y1))
# print((y2))
# ax = plt.subplot(111)
# ax.plot(ts, y2)
# ax2 = plt.twinx(ax)
# ax2.plot(ts, y1, color="r")

# ax.set_ylabel("Biomass", color="b")
# ax2.set_ylabel("Glucose", color="r")

# plt.show()
