luisa_genes = set( ['HP_0009', 'HP_0036', 'HP_0043', 'HP_0048', 'HP_0064', 'HP_0091', 'HP_0099', 'HP_0104', 'HP_0108', 'HP_0142', 'HP_0235', 'HP_0270', 'HP_0286', 'HP_0289', 'HP_0312', 'HP_0342', 'HP_0415', 'HP_0430', 'HP_0519', 'HP_0566', 'HP_0568', 'HP_0588', 'HP_0599', 'HP_0656', 'HP_0661', 'HP_0743', 'HP_0755', 'HP_0759', 'HP_0762', 'HP_0783', 'HP_0793', 'HP_0818', 'HP_0860', 'HP_0868', 'HP_0897', 'HP_0963', 'HP_0965', 'HP_1046', 'HP_1072', 'HP_1100', 'HP_1105', 'HP_1107', 'HP_1141', 'HP_1184', 'HP_1213', 'HP_1229', 'HP_1277', 'HP_1282', 'HP_1349', 'HP_1455', 'HP_1533'])

markus_genes = {'HP_0099', 'HP_0091', 'HP_1141', 'HP_0783', 'HP_1277', 'HP_0036', 'HP_0963', 'HP_0519', 'HP_0762', 'HP_0793', 'HP_0064', 'HP_0759', 'HP_0661', 'HP_1072', 'HP_0043', 'HP_0270', 'HP_1455', 'HP_1107', 'HP_0312', 'HP_0235', 'HP_0415', 'HP_1046', 'HP_1282', 'HP_0142', 'HP_1184', 'HP_0108', 'HP_0868', 'HP_1213', 'HP_1533', 'HP_0289', 'HP_0566', 'HP_0860', 'HP_0599', 'HP_0755', 'HP_0342', 'HP_0286', 'HP_0818', 'HP_1105', 'HP_0048', 'HP_0568', 'HP_0104', 'HP_0965', 'HP_0558', 'HP_0656', 'HP_0897', 'HP_1349', 'HP_1100', 'HP_0743', 'HP_0744_2', 'HP_0430', 'HP_1229'}


print(len(luisa_genes))
print(len(markus_genes))

print(len(markus_genes.intersection(luisa_genes)))

for x in sorted(luisa_genes.union(markus_genes)):

    if not x in luisa_genes or not x in markus_genes:
        print(x, x in luisa_genes, x in markus_genes)