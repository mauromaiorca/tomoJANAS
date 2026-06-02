## CryoSPARC integration

# setup
janas csparc_setup --license-id XXX  --host XXX --base-port 39000 --email XXX
source ~/.janas/cryosparc_env.sh



# setup selection project

--cs_project

# examples

janas_utils csparc_localnurefinement --particle-dir . --project P31 --workspace W1 --lane default --sym C1 --ref  "$(pwd)/class_1_rec.mrc" --mask mask_from_modelDilatedClose_softEdge.mrc --resplit class_1.star class_1_updatedLNU --precomputed J89
