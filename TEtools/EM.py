import numpy as np
import scipy.sparse as sparse


class expectation_maximization(object):
    def __init__(self, eff_lengths):
        self.min_step = 1.0
        self.max_step = 1.0
        self.inc_step = 4.0
        self.tolerance = 0.0001
        self.max_iterations = 100
        self.effective_lengths = eff_lengths

    def normalize(self, sparse_vect):
        tmp_sum = np.sum(sparse_vect.data)
        if tmp_sum > 0:
            norm = 1.0 / tmp_sum
        else:
            norm = 0.0
        return sparse_vect.multiply(norm)

    def compute_abundances(self, weight_vect, umi_masks):
        count_vect = sparse.csr_matrix(weight_vect.shape)
        for mask in umi_masks:
            tmp_cnt = mask.multiply(weight_vect)
            tmp_cnt = self.normalize(tmp_cnt)
            count_vect += tmp_cnt
        return count_vect

    def calc_step_weights(self, counts):
        weights = counts.multiply(self.effective_lengths)
        weights = self.normalize(weights)
        return weights

    def calc_diff_magnitude(self, new_weights, old_weights):
        difference = new_weights - old_weights
        magnitude = np.sqrt(np.sum((difference.multiply(difference)).data))
        return difference, magnitude

    def calc_alphaS(self, step_magnitude, course_correction):
        alphaS = step_magnitude / course_correction
        alphaS = max(self.min_step, min(self.max_step, alphaS))
        return alphaS

    def check_min_max(self, alphaS):
        if alphaS == self.max_step:
            self.max_step = self.inc_step * self.max_step

        if 0 > self.min_step == alphaS:  # check coverage, when would this condition occur?
            self.min_step = self.inc_step * self.min_step  # but we've shown that min_step is negative now?

    def preprocess(self, unmoored_counts, count_vect):
        umi_masks = []
        for count in unmoored_counts:
            count_vect += count
            umi_masks.append(count != 0)

        return sparse.csr_matrix(count_vect), umi_masks

    def run(self, unmoored_counts, gene_counts):
        count_vect, umi_masks = self.preprocess(unmoored_counts,
                                                np.zeros(max(self.effective_lengths.shape)))
        initial_weights = self.calc_step_weights(count_vect+gene_counts)

        iteration_cntr = 0

        while iteration_cntr < self.max_iterations:
            iteration_cntr += 1
            initial_counts = self.compute_abundances(initial_weights, umi_masks)
            first_step_weights = self.calc_step_weights(initial_counts+gene_counts)

            first_step_diff, first_step_magnitude = self.calc_diff_magnitude(first_step_weights, initial_weights)
            if first_step_magnitude < self.tolerance:
                break

            first_step_counts = self.compute_abundances(first_step_weights, umi_masks)
            second_step_weights = self.calc_step_weights(first_step_counts+gene_counts)

            second_step_diff, second_step_magnitude = self.calc_diff_magnitude(second_step_weights, first_step_weights)
            if second_step_magnitude < self.tolerance:
                initial_weights = second_step_weights
                break

            velocity, velocity_magnitude = self.calc_diff_magnitude(second_step_diff, first_step_diff)
            if velocity_magnitude < self.tolerance:
                initial_weights = first_step_weights  # confirm with Ying's math, seems off
                break

            course_correction = np.sqrt(abs(np.sum((first_step_diff.multiply(velocity)).data)))

            alphaS = self.calc_alphaS(first_step_magnitude, course_correction)

            new_weights = sparse.csr_matrix(
                np.fmax(0.0, initial_weights.todense()[0] + (2 * alphaS) * first_step_diff + (alphaS ** 2) * velocity))

            if abs(alphaS - 1.0) > 0.01:
                try:
                    new_counts = self.compute_abundances(new_weights, umi_masks)
                    new_weights = self.calc_step_weights(new_counts+gene_counts)
                except:
                    print("Error in EMupdate")
                    raise

            self.check_min_max(alphaS)

            new_weights, initial_weights = initial_weights, new_weights

        if iteration_cntr > self.max_iterations:
            print("did not converge by " + str(self.max_iterations) + " iterations. \n")

        return self.compute_abundances(initial_weights, umi_masks)


def EM(unmo_vect, effective_lengths, gene_counts):
    em_loops = expectation_maximization(effective_lengths)
    unmo_cnt = em_loops.run(unmo_vect, gene_counts)
    return unmo_cnt

