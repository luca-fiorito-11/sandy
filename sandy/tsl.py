# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for Thermal Neutron
Scattering Data.

class `Tsl` class acts as a container for temperatura-dependent tabulated
thermal neutron scattering cross section.
"""

import pandas as pd
import sandy
import pytest 

class Tsl():
    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, dct):
        self.data = dct

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = dict(data)

    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract hierarchical structure of tsl data from `sandy.Endf6`
        instance.

        Parameters
        ----------
        tape : `sandy.Endf6`
            instance containing tsl data

        Returns
        -------
        `sandy.Tsl`
            structured container divided into the 3 reaction of the tsl files:
            elastic coherent scattering, elastic incoherent scattering and
            inelastic incoherent scattering.

        Examples
        --------
        Load test ENDF-6 file with data for Be-4:
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 26)
        >>> tsl = Tsl.from_endf6(tape)

        Coherent elastic scattering:
        >>> tsl.data['elastic coherent']['T'].keys()
        dict_keys([296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0])

        Incoherent inelastic scattering:
        >>> tsl.data['inelastic incoherent']['beta'][0.0]['T'].keys()
        dict_keys([296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0])

        Incoherent elastic scattering:
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 10)
        >>> tsl = Tsl.from_endf6(tape)
        >>> tsl.data['elastic incoherent']['Debye-Waller']
        array([14.70372, 19.1224 , 20.37892, 21.65261, 21.97355, 22.94205,
               23.26671, 24.24591, 24.57398])
        """
        tape = endf6.filter_by(listmf=[7], listmt=[2, 4])
        data = {}
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            if mt == 2:
                LTHR = sec['LTHR']
                if LTHR == 1:  # Coherent elastic scattering
                    add = {'E': sec["EINT"]}
                    add_2 = {}
                    for temp, S in sec['T'].items():
                        add_2[temp] = {'S': S['S']}
                    add['T'] = add_2
                    data['elastic coherent'] = add
                elif LTHR == 2:  # Incoherent elastic scattering
                    add = {
                        'characteristic bound xs': sec['SB'],
                        'T': sec['TINT'],
                        'Debye-Waller': sec['W'],
                    }
                    data['elastic incoherent'] = add
            elif mt == 4:  # Inelastic incoherent scattering
                add = {
                    'alpha': sec["alpha"]
                }
                add_2 = {}
                for beta, temp_info in sec['beta/T'].items():
                    add_3 = {}
                    for temp, S_info in temp_info.items():
                        add_3[temp] = {'S': S_info['S']}
                    add_2[beta] = {'T': add_3}
                add['beta'] = add_2
                data['inelastic incoherent'] = add
        return cls(data)

    def _S_elastic_coherent(self):
        """
        The S matrix for elastic coherent scattering.

        Returns
        -------
        S : `list`
            List of dataframes containing the stack S matrix.

        Examples
        --------
        Load test ENDF-6 file with data for Be-4:
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 26)
        >>> tsl = Tsl.from_endf6(tape)
        >>> tsl._S_elastic_coherent()[0][0:5]
                      T	          E	          S
        0	2.96000e+02	1.62650e-03	0.00000e+00
        1	2.96000e+02	5.28450e-03	8.69813e-03
        2	2.96000e+02	6.50600e-03	1.91051e-02
        3	2.96000e+02	6.91100e-03	6.40282e-02
        4	2.96000e+02	1.17905e-02	7.49635e-02
        """
        data = self.data['elastic coherent']
        S = []
        E = data['E']
        for temp, s_info in data['T'].items():
            S_part = s_info['S']
            index = pd.MultiIndex.from_product([[temp], E],
                                               names=['T', 'E'])
            df = pd.DataFrame(S_part, index=index, columns=['S']).stack()\
                   .reset_index().drop(['level_2'], axis=1)\
                   .rename(columns={0: 'S'})
            S.append(df)
        return S

    def _S_inelastic_incoherent(self):
        """
        The S matrix for inelastic incoherent scattering.

        Returns
        -------
        S : `list`
            List of dataframes containing the stack S matrix.

        Examples
        --------
        Load test ENDF-6 file with data for Be-4:
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 26)
        >>> tsl = Tsl.from_endf6(tape)
        >>> tsl._S_inelastic_incoherent()[0][0:5]
                   beta	          T	      alpha	          S
        0	0.00000e+00	2.96000e+02	3.05297e-03	7.52844e-05
        1	0.00000e+00	2.96000e+02	3.27092e-03	8.06581e-05
        2	0.00000e+00	2.96000e+02	3.50442e-03	8.64153e-05
        3	0.00000e+00	2.96000e+02	3.75460e-03	9.25834e-05
        4	0.00000e+00	2.96000e+02	4.02263e-03	9.91918e-05
        """
        data = self.data['inelastic incoherent']
        S = []
        alpha = data['alpha']
        for beta, t_info in data['beta'].items():
            for temp, s_info in t_info['T'].items():
                S_part = s_info['S']
                index = pd.MultiIndex.from_product([[beta], [temp], alpha],
                                                   names=['beta', 'T',
                                                          'alpha'])
                df = pd.DataFrame(S_part, index=index, columns=['S'])\
                       .stack().reset_index().drop(['level_3'], axis=1)\
                       .rename(columns={0: 'S'})
                S.append(df)
        return S

    def get_S(self, kind='elastic coherent'):
        """
        Extract stacked S matrix from the Tsl class.

        Parameters
        ----------
        kind : `str`, optional
            The type of reaction from where the S matrix is going to be
            extracted. The default is 'elastic coherent'.

        Notes
        -----
        .. note:: In a future, the return of this method is going to be a
                  object of the class `S_matrix`.

        Returns
        -------
        `pd.DataFrame`
            Stack S matrix.

        Examples
        --------
        Load test ENDF-6 file with data for Be-4:
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 26)
        >>> tsl = Tsl.from_endf6(tape)

        Elastic coherent S-matrix:
        >>> tsl.get_S(kind ='elastic coherent').head()
                      T	          E	          S
        0	2.96000e+02	1.62650e-03	0.00000e+00
        1	2.96000e+02	5.28450e-03	8.69813e-03
        2	2.96000e+02	6.50600e-03	1.91051e-02
        3	2.96000e+02	6.91100e-03	6.40282e-02
        4	2.96000e+02	1.17905e-02	7.49635e-02

        Inelastic incoherent S-matrix:
        >>> tsl.get_S(kind ='inelastic incoherent').head()
                   beta	          T	      alpha	          S
        0	0.00000e+00	2.96000e+02	3.05297e-03	7.52844e-05
        1	0.00000e+00	2.96000e+02	3.27092e-03	8.06581e-05
        2	0.00000e+00	2.96000e+02	3.50442e-03	8.64153e-05
        3	0.00000e+00	2.96000e+02	3.75460e-03	9.25834e-05
        4	0.00000e+00	2.96000e+02	4.02263e-03	9.91918e-05
        """
        if kind == 'elastic coherent':
            S = self._S_elastic_coherent()
        elif kind == 'inelastic incoherent':
            S = self._S_inelastic_incoherent()
        return pd.concat(S)

    def to_endf6(self, endf6):
        """
        Update S matrixs in `Endf6` instance with those available in a
        `Tsl` instance.

        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance

        Notes
        -----
        .. note:: The user can only add new energy grids in the coherent
                  elastic scattering. For the moment, add new temperatures is
                  not allow.
        .. note:: The user can only add new alpha grids in the incoherent
                  inelastic scattering. For the moment, add new temperatures
                  and beta values is not allow.

        Returns
        -------
        `sandy.Endf6`
            `Endf6` instance with updated S matrixs and Deby-Waller function

        Examples
        --------
        Incoherent elastic scattering
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 10)
        >>> from_endf = sandy.sections.mf7.read_mf7(tape, 10, 2)
        >>> text = sandy.sections.mf7.write_mf7(from_endf)
        >>> tsl = sandy.Tsl.from_endf6(tape)
        >>> new_tape = tsl.to_endf6(tape)
        >>> new_from_endf = sandy.sections.mf7.read_mf7(tape, 10, 2)
        >>> new_text = sandy.sections.mf7.write_mf7(new_from_endf)
        >>> assert new_text == text

        Coherent elastic scattering
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 26)
        >>> from_endf = sandy.sections.mf7.read_mf7(tape, 26, 2)
        >>> text = sandy.sections.mf7.write_mf7(from_endf)
        >>> tsl = sandy.Tsl.from_endf6(tape)
        >>> new_tape = tsl.to_endf6(tape)
        >>> new_from_endf = sandy.sections.mf7.read_mf7(tape, 26, 2)
        >>> new_text = sandy.sections.mf7.write_mf7(new_from_endf)
        >>> assert new_text == text

        Incoherent inelastic scattering
        >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', 26)
        >>> from_endf = sandy.sections.mf7.read_mf7(tape, 26, 4)
        >>> text = sandy.sections.mf7.write_mf7(from_endf)
        >>> tsl = sandy.Tsl.from_endf6(tape)
        >>> new_tape = tsl.to_endf6(tape)
        >>> new_from_endf = sandy.sections.mf7.read_mf7(tape, 26, 4)
        >>> new_text = sandy.sections.mf7.write_mf7(new_from_endf)
        >>> assert new_text == text
        """
        data_endf6 = endf6.data.copy()
        tape = endf6.filter_by(listmf=[7], listmt=[2, 4])
        for (mat, mf, mt) in tape.keys:
            sec = tape.read_section(mat, mf, mt)
            if mt == 2:
                LTHR = sec['LTHR']
                if LTHR == 1:
                    obj_data = self.data['elastic coherent']
                    sec["EINT"] = obj_data['E']
                    for temp, S_info in sec['T'].items():
                        sec['T'][temp]['S'] = obj_data['T'][temp]['S']
                elif LTHR == 2:
                    obj_data = self.data['elastic incoherent']
                    sec['SB'] = obj_data['characteristic bound xs']
                    sec['TINT'] = obj_data['T']
                    sec['W'] = obj_data['Debye-Waller']
            elif mt == 4:
                obj_data = self.data['inelastic incoherent']
                sec['alpha'] = obj_data['alpha']
                for beta, beta_info in sec['beta/T'].items():
                    obj_beta = obj_data['beta'][beta]
                    for temp, temp_info in beta_info.items():
                        sec['beta/T'][beta][temp]['S'] = obj_beta['T'][temp]['S']
            data_endf6[(mat, mf, mt)] = sandy.write_mf7(sec)
        return sandy.Endf6(data_endf6)
