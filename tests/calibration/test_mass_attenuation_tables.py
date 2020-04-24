'''Test the mass attenuation tables'''

import unittest
from ogo.calibration.mass_attenuation_tables import \
    mass_attenuation_tables


class TestMassAttenuationTables(unittest.TestCase):
    '''Test the mass attenuation tables'''

    def runner(self, name, length):
        '''Run a test'''
        self.assertTrue(name in mass_attenuation_tables)
        table = mass_attenuation_tables[name]

        self.assertTrue(hasattr(table, 'Energy [keV]'))
        self.assertTrue(hasattr(table, 'Mass Attenuation [cm2/g]'))

        self.assertEqual(table.shape[0], length)

    def test_adipose_table(self):
        '''Adipose table is as expected'''
        self.runner('adipose_table', 44)

    def test_air_table(self):
        '''Air table is as expected'''
        self.runner('air_table', 38)

    def test_blood_table(self):
        '''Blood table is as expected'''
        self.runner('blood_table', 51)

    def test_bone_table(self):
        '''Bone table is as expected'''
        self.runner('bone_table', 49)

    def test_muscle_table(self):
        '''Muscle table is as expected'''
        self.runner('muscle_table', 49)

    def test_k2hpo4_table(self):
        '''K2HPO4 table is as expected'''
        self.runner('k2hpo4_table', 50)

    def test_cha_table(self):
        '''CHA table is as expected'''
        self.runner('cha_table', 50)

    def test_triglyceride_table(self):
        '''Triglyceride table is as expected'''
        self.runner('triglyceride_table', 46)

    def test_water_table(self):
        '''Water table is as expected'''
        self.runner('water_table', 36)


if __name__ == '__main__':
    unittest.main()
