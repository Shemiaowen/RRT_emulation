# 90d risk difference
data1 <- data.frame(
  Trial = c("Target trial 1", "Target trial 2"),
  RiskDifference = c(0.026, -0.017),
  LowerCI = c(0.022, -0.025),
  UpperCI = c(0.03, -0.008)
)

# forest plot
ggplot(data1, aes(x = Trial, y = RiskDifference)) +
  geom_point(color = "red", size = 2) +  
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05, color = "black") +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  
  labs(y = "90-d risk difference",  
       x = "") +  
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  
  scale_y_continuous(breaks = seq(-0.05, 0.05, by = 0.01))

# 30d risk difference
data2 <- data.frame(
  Trial = c("Target trial 1", "Target trial 2"),
  RiskDifference = c(0.018, -0.025),
  LowerCI = c(0.014, -0.033),
  UpperCI = c(0.022, -0.016)
)

# forest plot
ggplot(data2, aes(x = Trial, y = RiskDifference)) +
  geom_point(color = "red", size = 2) +  
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05, color = "black") +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  
  labs(y = "30-d risk difference",  
       x = "") +  
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  
  scale_y_continuous(breaks = seq(-0.05, 0.05, by = 0.01))
