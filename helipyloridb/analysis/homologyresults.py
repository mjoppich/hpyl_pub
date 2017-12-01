from database.homologydb import HomologyDatabase


class HomologyResult:

    def __init__(self):

        self.homology_relations = []
        self.combination_results = []
        self.mul_combination_results = []

    def toDataBase(self, homolDB):

        if self.homology_relations != None:

            for relation in self.homology_relations:
                homolDB.addHomologyRelation( relation[0], relation[1], relation[2] )

        if self.combination_results != None:

            for comb_result in self.combination_results:
                homolDB.addCombination( comb_result[0], comb_result[1], comb_result[2] )

        if self.mul_combination_results != None:

            for mul_comb_result in self.mul_combination_results:
                homolDB.addMultiCombination(mul_comb_result)