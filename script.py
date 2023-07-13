import cbmpy
from cbmpy.CBModel import Model
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.DynamicModels.DynamicJointFBA import DynamicJointFBA

m_a: Model = model_a.build_joint_fba_model_A()
m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_3").setUpperBound(3)
m_a.getReaction("R_5").setUpperBound(1)


m_b: Model = model_b.build_joint_fba_model_B()
m_b.getReaction("R_1").setUpperBound(10)
m_b.getReaction("R_3").setUpperBound(1)
m_b.getReaction("R_5").setUpperBound(1)

community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)

dj = DynamicJointFBA(community_model, [1, 1])
T, metabolites, biomasses = dj.simulate(0.1)

print(T)
print(metabolites)
print(biomasses)
