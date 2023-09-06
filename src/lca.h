inline float
updateLCAtProbabilityC2N(std::vector<float>& p_vec)
{
  float prob1 = 1;
  float prob2 = 0;
  for (unsigned int i = 0; i < p_vec.size(); ++i)
    prob1 = prob1 * (1 - p_vec[i]);
  for (unsigned int i = 0; i < p_vec.size(); ++i) {
    float p_i = p_vec[i];
    for (unsigned int j = 0; j < p_vec.size(); ++j) {
      if (j != i)
        p_i = p_i * (1 - p_vec[j]);
    }
    prob2 += p_i;
  }
  return 1 - prob1 - prob2;
}

inline float
updateLCAtProbabilityAPN(std::vector<float>& p_vec, std::vector<float>& q_vec)
{
  assert(p_vec.size() == q_vec.size());
  float p_term = 0;
  for (unsigned int i = 0; i < p_vec.size(); ++i) {
    float p_in = 0;
    for (unsigned int j = 0; j < p_vec.size(); ++j) {
      if (i != j) {
        p_in += p_vec[i] * p_vec[j] * q_vec[j] / (1 - q_vec[i]);
      }
    }
    p_term += p_in * q_vec[i];
  }
  return 1 - pow(1 - p_term, p_vec.size() * q_vec.size());
}
