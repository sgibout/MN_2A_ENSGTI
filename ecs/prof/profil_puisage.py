import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def generate_puisage_profile(num_people=4, time_step_minutes=6):
    # Durée totale : 24 heures, pas de temps : 6 minutes
    time_steps = 24 * 60 // time_step_minutes
    total_seconds = 24 * 60 * 60
    time_step_seconds = time_step_minutes * 60

    # Initialisation du profil de puisage (en litres)
    profile = np.zeros(time_steps)

    # Matin (6h-9h) : Douches, lavage des mains, petit-déjeuner
    profile[60:90] = np.random.normal(30, 10, 30) * num_people / 4

    # Après-midi (12h-14h) : Cuisine, lavage des mains
    profile[120:140] = np.random.normal(10, 3, 20) * num_people / 4

    # Soir (18h-22h) : Douches, dîner, lavage de la vaisselle
    profile[180:220] = np.random.normal(25, 8, 40) * num_people / 4

    # Assurer que les valeurs restent positives
    profile = np.clip(profile, 0, None)

    # Génération du temps en secondes à partir de 0
    time_seconds = np.arange(0, total_seconds, time_step_seconds)

    # Création du DataFrame pour une meilleure visualisation
    df = pd.DataFrame({"Temps (s)": time_seconds, "Consommation (L)": profile})

    # Création du DataFrame pour une meilleure visualisation
    #time_index = pd.date_range(start="00:00", end="23:54", freq=f"{time_step_minutes}min")
    #df = pd.DataFrame(profile, index=time_index, columns=["Consommation (L)"])

    # Sauvegarde du profil sous forme de fichier CSV

    return df

num_people=4
time_step_minutes=6
df = generate_puisage_profile(num_people=num_people,time_step_minutes=time_step_minutes)

csv_path = 'profil_puisage.csv'
df.to_csv(csv_path)

# Affichage du profil
plt.figure(figsize=(12, 6))
plt.plot(df.index, df["Consommation (L)"], label="Profil de puisage")
plt.title(f"Profil de puisage d'eau chaude pour une famille de {num_people} personnes (Pas de {time_step_minutes} minutes)")
plt.xlabel("Temps")
plt.ylabel("Consommation (L)")
plt.grid(True)
plt.xticks(rotation=45)
plt.legend()
plt.show()

# Appel de la fonction pour une famille de 4 personnes
