%% script to show the lookuptable of MP2RAGE


struct_MP2RAGE.ETL = 128;
struct_MP2RAGE.TI1 = 800;
struct_MP2RAGE.TI2 = 2000;
struct_MP2RAGE.alpha1 = 4;
struct_MP2RAGE.alpha2 = 7;
struct_MP2RAGE.MP2RAGE_TR = 5000;
struct_MP2RAGE.TR = 7;

[struct_out] = MP2RAGE_LookUpTable_Only(struct_MP2RAGE)