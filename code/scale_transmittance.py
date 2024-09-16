def scale_transmittance(final_transmittance, Q_target_old, Q_target_new):
    scaling_factor = Q_target_new / Q_target_old
    
    # Scale the transmittance
    final_transmittance_new = final_transmittance ** scaling_factor
    
    return final_transmittance_new