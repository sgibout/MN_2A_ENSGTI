import numpy as np


class State:
    def __init__(self):
        self.T = np.nan  # Température [°C]
        self.P = np.nan  # Pression [bar]
        self.u = np.nan  # énergie interne massique [J/kg]
        self.h = np.nan  # enthalpie massique [J/kg]
        self.s = np.nan  # entropie massique [J/K/kg]
        self.v = np.nan  # volume spécifique [m3/kg]

    def set_T_P(self, T, P):
        """
        Calcule l'état complet en fonction de la température et de la pression
        :param P: Pression
        :param T: Température
        :return: l'état complet
        """
        self.P = P
        self.T = T
        self.u = State.u_T_P(T, P)
        self.h = State.h_T_P(T, P)
        self.s = State.s_T_P(T, P)
        self.v = State.v_T_P(T, P)

        return self

    def set_h_P(self, h, P):
        """
        Calcule l'état complet en fonction de l'enthalpie massique et de la pression
        :param h: enthalpie massique
        :param T: Température
        :return: l'état complet
        """
        self.P = P
        self.T = T = State.T_h_P(h, P)
        self.u = State.u_T_P(T, P)
        self.h = h
        self.s = State.s_T_P(T, P)
        self.v = State.v_T_P(T, P)

        return self

    def set_u_v(self, u, v):
        """
        Calcule l'état complet en fonction de l'énergie (interne) massique et du volume massique
        :param u: enthalpie massique
        :param v: Température
        :return: l'état complet
        """
        self.u = u
        self.T = T = State.T_u_P(u, np.nan) # pas besoin de la pression ici

        self.P = P = State.P_T_v(T,v)
        self.h = State.h_T_P(T, P)
        self.s = State.s_T_P(T, P)
        self.v = State.v_T_P(T, P)

        return self

    def copy(self, state):
        self.T = state.T
        self.P = state.P
        self.u = state.u
        self.h = state.h
        self.s = state.s
        self.v = state.v

        return self

    # les méthodes statiques (et constantes)
    # Etat de référence
    T0 = 273.15  # K
    P0 = 1E5  # Pa
    u0 = 0
    h0 = 0
    s0 = 0
    # Constantes
    CP = 1000
    CV = 713.31
    R = 8.314  # J/mol/K
    Mmol = 29E-3  # kg/mol
    r = R/Mmol # J/kg/K

    def __str__(self):
        return f"Etat T={self.T} P={self.P} u={self.u} h={self.h} s={self.s} v={self.v}"

    @staticmethod
    def h_T_P(T, P):
        h = State.h0 + State.CP * (T + 273.15 - State.T0)
        return h

    @staticmethod
    def T_h_P(h, P):
        T = State.T0 - 273.15 + (h - State.h0) / State.CP
        return T

    @staticmethod
    def u_T_P(T, P):
        u = State.u0 + State.CV * (T + 273.15 - State.T0)
        return u

    @staticmethod
    def T_u_P(u, P):
        T = State.T0 - 273.15 + (u - State.u0) / State.CV
        return T

    @staticmethod
    def s_T_P(T, P):
        s = State.s0 + State.CP * np.log((T + 273.15) / State.T0) + State.CV * np.log(1e5 * P / State.P0)
        return s

    @staticmethod
    def v_T_P(T, P):
        v = State.R * (T + 273.15) / State.Mmol / (1E5 * P)
        return v

    @staticmethod
    def P_T_v(T, v):
        P = State.R * (T + 273.15) / (1E5 * v * State.Mmol)
        return P
