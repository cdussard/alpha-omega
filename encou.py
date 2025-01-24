# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:26:50 2024

@author: claire.dussard
"""
    
    ind = np.where(freq_med==55)[0][0]

    ind_beta = np.where(freq_med==35)[0][0]
    
    psd_corr = psd_med / np.sum(psd_med[ind:])
    freq_med_ = freq_med[0:ind_beta]
    psd_corr = psd_corr[0:ind_beta]
    ax[1].plot(freq_med_, psd_corr, label=descr)
    
    
    
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=False)

# Indices pour les régions "Central" et "Posterior"
el_c = np.where([desc_[0] == 'Central' for desc_ in desc])[0]
files_c = [files[i] for i in el_c]
desc_c = [desc[i] for i in el_c]

el_p = np.where([desc_[0] == 'Posterior' for desc_ in desc])[0]
files_p = [files[i] for i in el_p]
desc_p = [desc[i] for i in el_p]

# Boucle pour les signaux "Central" (gauche)
for i, (sig, descr) in enumerate(zip(files_c, desc_c)):
    print(sig)
    print(descr)
    times = create_times(len(sig) / fs, fs)

    # Vérifier et ajuster la longueur du signal
    if len(times) != len(sig):
        min_length = min(len(times), len(sig))
        times = times[:min_length]
        sig = sig[:min_length]

    # Calcul du spectre
    freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=int(fs * 1.2), f_range=[fmin, fmax])
    ind_beta = np.where(freq_med == 35)[0][0]

    # Traçage sur le sous-plot de gauche
    axes[0].plot(freq_med[:ind_beta], psd_med[:ind_beta], label=descr)

# Boucle pour les signaux "Posterior" (droite)
for i, (sig, descr) in enumerate(zip(files_p, desc_p)):
    print(sig)
    print(descr)
    times = create_times(len(sig) / fs, fs)

    # Vérifier et ajuster la longueur du signal
    if len(times) != len(sig):
        min_length = min(len(times), len(sig))
        times = times[:min_length]
        sig = sig[:min_length]

    # Calcul du spectre
    freq_med, psd_med = compute_spectrum(sig, fs, method='welch', avg_type='median', nperseg=int(fs * 1.2), f_range=[fmin, fmax])
    ind_beta = np.where(freq_med == 35)[0][0]

    # Traçage sur le sous-plot de droite
    axes[1].plot(freq_med[:ind_beta], psd_med[:ind_beta], label=descr)

# Ajout des légendes
axes[0].legend(title="Central (Uncorrected)")
axes[1].legend(title="Posterior (Corrected)")

# Ajustement de la mise en page
plt.tight_layout()
plt.show()
