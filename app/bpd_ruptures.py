import ruptures as rpt

def analyze_ruptures_breakpoints(original_signal, start, end, penalty=15000, min_size=4):
    # Slice the signal according to the start and end indices
    signal_segment = original_signal[start:end]
    
    # Initialize the PELT algorithm with the l2 model and minimum segment size
    algo = rpt.Pelt(model="l2", min_size=min_size)
    
    # Fit the model to the signal segment
    algo.fit(signal_segment)
    
    # Predict the breakpoints with the given penalty
    result = algo.predict(pen=penalty)
    
    # Adjust the breakpoints to reflect the position in the original signal
    adjusted_result = [r + start for r in result]
    
    # Calculate the number of breakpoints
    number_of_breakpoints = len(adjusted_result) - 1
    
    # # Display the signal with the breakpoints
    # plt.figure(figsize=(9.15, 4))
    # rpt.display(signal_segment, [], result)
    # plt.show()
    
    return adjusted_result