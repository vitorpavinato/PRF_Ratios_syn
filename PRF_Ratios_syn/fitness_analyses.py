"""
Module containing function to estimate fitness for amino acids with two or more codons.
"""

import itertools
from typing import Callable, Any, Tuple, List, Dict
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
from scipy import stats
from sklearn.linear_model import LinearRegression


def extract_unique_codons(codon_changes: List[Tuple[str, float]]) -> List[str]:
    """
    Extract unique codons from a list of codon changes.
    """
    unique_codons = set()
    for change, _ in codon_changes:
        start, end = change.split('>')
        unique_codons.add(start)
        unique_codons.add(end)

    return sorted(list(unique_codons))


def import_from_codonpair(inputfile: str) -> Dict[str, Tuple[List[Tuple[str, float]], List[str]]]:
    """
    Import codon pair data from a file and extract unique codons for each amino acid.
    """
    
    codon_data: Dict[str, Tuple[List[Tuple[str, float]], List[str]]] = {}

    try:
        with open(inputfile, 'r', encoding='utf-8') as infile:
            for line in infile:
                if line.startswith("Codon_pair"):
                    continue
                if len(line.split()) == 8 and ">" in line.split()[0]:
                    parts = line.split()
                    codon_change = parts[0]
                    amino_acid = parts[1]
                    try:
                        two_ns = float(parts[7])
                    except ValueError:
                        continue  # Skip lines where conversion to float fails

                    # Add the codon change to the dictionary
                    if amino_acid not in codon_data:
                        codon_data[amino_acid] = ([],[])

                    codon_data[amino_acid][0].append((codon_change, two_ns))
        
        # After processing all lines, calculate unique codons for each amino acid
        for amino_acid, (changes, _) in codon_data.items():
            print(changes)
            unique_codons = extract_unique_codons(changes)
            codon_data[amino_acid] = (changes, unique_codons)

    except FileNotFoundError:
        print(f"Error: File '{inputfile}' not found.")
    except IOError:
        print(f"Error: Could not read file '{inputfile}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return codon_data


def create_codon_matrix(
        codon_changes: List[Tuple[str, float]],
        codons: List[str]
        ) -> np.ndarray[np.float64]:
    """
    Function to make n x n matrix of codon changes, being the n the number of codons.
    The n x n matrix will be compatible with downstream linear regression model.
    """

    # Create a 4x4 matrix of zeros
    m = np.zeros((len(codons), len(codons)))

    # Create a dictionary to map codon pairs to their 2Ns values
    codon_dict = {pair: value for pair, value in codon_changes}

    # Populate the matrix
    for i, from_codon in enumerate(codons):
        for j, to_codon in enumerate(codons):
            if i != j:  # Skip diagonal elements
                pair = f"{from_codon}>{to_codon}"
                reverse_pair = f"{to_codon}>{from_codon}"

                if pair in codon_dict:
                    m[i][j] = codon_dict[pair]
                elif reverse_pair in codon_dict:
                    m[i][j] = -codon_dict[reverse_pair]
                # If neither pair is in the dictionary, it remains 0

    # Return the matrix
    return m


def build_design_matrix(num_codons: int) -> np.ndarray:
    """
    Build the design matrix A for the linear model Ax = b.
    
    :param num_codons: Number of codons (minimum 3)
    :return: Design matrix A
    """

    # Check that the number of codons is at least 3
    if num_codons < 3:
        raise ValueError("Number of codons must be at least 3")

    # Generate all pairs of codon indices
    pairs = list(itertools.combinations(range(num_codons), 2))

    # Initialize the design matrix
    A = np.zeros((len(pairs), num_codons - 1))

    for i, (j, k) in enumerate(pairs):
        if j == 0:  # First codon is reference (s_1 = 0)
            A[i, k-1] = -1
        elif k == 0:  # First codon is reference (s_1 = 0)
            A[i, j-1] = 1
        else:
            A[i, j-1] = 1
            A[i, k-1] = -1

    return A


def build_observation_vector(m: np.ndarray) -> np.ndarray:
    """
    Build the observation vector b for the linear model Ax = b,
    using only the upper triangular part of the matrix.
    
    :param m: The matrix of codon changes
    :return: Observation vector b
    """
    
    num_codons = m.shape[0]
    pairs = list(itertools.combinations(range(num_codons), 2))

    # Extract only upper triangular values
    b = np.array([m[i, j] for i, j in pairs])

    return b


def solve_least_squares(A: np.ndarray, b: np.ndarray) -> dict:
    """
    Solve the least squares problem for fitness values and calculate various metrics.
    
    :param A: Design matrix
    :param b: Observation vector
    :return: Dictionary containing various results and metrics
    """
    
    # Solve the least squares problem
    # x: the solution vector - the selection coefficients
    # residual_ss: the residual sum of squares = sum((b - Ax) ** 2)
    # rank: the rank of the matrix
    # s: the singular values
    x, residual_ss, rank, s = lstsq(A, b)

    # Extract the results
    p_vector = np.concatenate(([1], 1 + x))

    # Calculate predicted values
    b_pred = A.dot(x)

    # Calculate error metrics
    relative_error_square = np.sum(((b - b_pred) / b_pred) ** 2)
    normalized_squared_error = np.sum(((b - b_pred) ** 2) / b_pred)

    # Calculate correlation coefficient
    # correlation_coefficient = np.corrcoef(b, b_pred)[0, 1]
    # Or if you prefer stats.pearsonr:
    correlation_coefficient, p_value = stats.pearsonr(b, b_pred)

    # Calculate R-squared and Adjusted R-squared
    total_ss = np.sum((b - np.mean(b)) ** 2)
    r_squared = 1 - (residual_ss / total_ss)
    adjusted_r_squared = 1 - (1 - r_squared) * (len(b) - 1) / (len(b) - len(x) - 1)

    # Calculate MSE, AIC and BIC
    mse = np.mean((b - b_pred) ** 2)
    rmse = np.sqrt(mse)
    aic = 2 * len(x) + len(b) * np.log(mse)
    bic = np.log(len(b)) * len(x) + len(b) * np.log(mse)

    # Calculate confidence intervals
    dof = max(0, A.shape[0] - A.shape[1])
    t_value = stats.t.ppf(0.975, dof)
    se = np.sqrt(np.sum((b - b_pred) ** 2) / dof * np.diag(np.linalg.inv(A.T.dot(A))))
    ci_lower = x - t_value * se
    ci_upper = x + t_value * se

    return {
        "rank": rank,
        "singular_values": s,
        "relative_fitness": p_vector,
        "observed_selection_coefficients": b,
        "predicted_selection_coefficients": b_pred,
        "residual_sum_of_squares": residual_ss,
        "relative_error_square": relative_error_square,
        "normalized_squared_error": normalized_squared_error,
        "correlation_coefficient": correlation_coefficient,
        "r_squared": r_squared,
        "adjusted_r_squared": adjusted_r_squared,
        "mse": mse,
        "rmse": rmse,
        "aic": aic,
        "bic": bic,
        "confidence_intervals": list(zip(ci_lower, ci_upper)),
        "degrees_of_freedom": dof
    }


# If plot without lower matrix values - don't need this
# But keep it here for now.
def construct_full_predicted_matrix(b_pred: np.ndarray) -> np.ndarray:
    """
    Construct a full 4x4 matrix of predicted b values from the b_pred vector.
    
    :param b_pred: Vector of predicted b values from A.dot(x)
    :return: Full 4x4 matrix of predicted b values
    """
    
    # Create a 4x4 matrix of zeros
    m_pred = np.zeros((4, 4))

    # Fill the upper triangle of the matrix
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            m_pred[i, j] = b_pred[idx]
            m_pred[j, i] = -b_pred[idx]  # Fill lower triangle with negative values
            idx += 1

    return m_pred


# m matrices version
def scatterplot_selcoef_matrices(
        codons: List[str],
        observed_selcoefs: np.ndarray,
        predicted_selcoefs: np.ndarray
        ) -> None:
    """
    Function to plot the observed and predicted selection coefficients.
    You should supply m and m_predicted as matrices of the same dimensions.

    Parameters:
    codons (list): List of codons
    observed_selcoefs (array-like matrix): Observed selection coefficients
    predicted_selcoefs (array-like matrix): Predicted selection coefficients
    """

    # Flatten and remove diagonal zeros
    observed_selcoefs_flat = observed_selcoefs[observed_selcoefs != 0]
    predicted_selcoefs_flat = predicted_selcoefs[predicted_selcoefs != 0]

    # Create the scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(observed_selcoefs_flat, predicted_selcoefs_flat, alpha=0.7)
    plt.xlabel('Observed Selection Coefficients')
    plt.ylabel('Predicted Selection Coefficients')
    plt.title('Observed vs. Predicted Selection Coefficients')

    # Add the perfect fit line (y=x)
    min_val = min(plt.xlim()[0], plt.ylim()[0])
    max_val = max(plt.xlim()[1], plt.ylim()[1])
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='y=x')

    # Fit a regression line
    reg = LinearRegression().fit(observed_selcoefs_flat.reshape(-1, 1), predicted_selcoefs_flat)
    plt.plot(observed_selcoefs_flat, reg.predict(observed_selcoefs_flat.reshape(-1, 1)), 'g-', label='Fitted Line')

    # Add labels for each point
    for i, codon_i in enumerate(codons):
        for j, codon_j in enumerate(codons):
            if i != j:
                # print(f'{codon_i}>{codon_j}')
                plt.annotate(f'{codon_i}>{codon_j}',
                            (observed_selcoefs[i][j], predicted_selcoefs[i][j]),
                             xytext=(5, 5), textcoords='offset points')


    plt.legend()
    plt.grid(False)

    # Add text with R-squared value
    r_squared = reg.score(observed_selcoefs_flat.reshape(-1, 1), predicted_selcoefs_flat)
    plt.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes, 
            verticalalignment='top')

    plt.tight_layout()
    plt.show()


# b vectors version
def scatterplot_selcoef_vectors(
        codons: List[str],
        observed_selcoefs: np.ndarray,
        predicted_selcoefs: np.ndarray
        ) -> None:
    """
    Function to plot the observed and predicted selection coefficients.
    You should supply b and b_predicted as vectors of the same length.

    Parameters:
    codons (list): List of codons
    observed_selcoefs (array-like vector): Observed selection coefficients
    predicted_selcoefs (array-like vector): Predicted selection coefficients
    """

    # Flatten and remove diagonal zeros
    observed_selcoefs_flat = observed_selcoefs[observed_selcoefs != 0]
    predicted_selcoefs_flat = predicted_selcoefs[predicted_selcoefs != 0]

    # Create the scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(observed_selcoefs_flat, predicted_selcoefs_flat, alpha=0.7)
    plt.xlabel('Observed Selection Coefficients')
    plt.ylabel('Predicted Selection Coefficients')
    plt.title('Observed vs. Predicted Selection Coefficients')

    # Add the perfect fit line (y=x)
    min_val = min(plt.xlim()[0], plt.ylim()[0])
    max_val = max(plt.xlim()[1], plt.ylim()[1])
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='y=x')

    # Fit a regression line
    reg = LinearRegression().fit(observed_selcoefs_flat.reshape(-1, 1), predicted_selcoefs_flat)
    plt.plot(observed_selcoefs_flat, reg.predict(observed_selcoefs_flat.reshape(-1, 1)), 'g-', label='Fitted Line')

    # Add labels for each point
    pairs = list(itertools.combinations(range(len(codons)), 2))
    codon_pairs = np.array([f'{codons[i]}>{codons[j]}' for i, j in pairs])
    
    for i, codon_pair in enumerate(codon_pairs):
        plt.annotate(codon_pair,
                    (observed_selcoefs_flat[i], predicted_selcoefs_flat[i]),
                     xytext=(5, 5), textcoords='offset points')

    plt.legend()
    plt.grid(False)

    # Add text with R-squared value
    r_squared = reg.score(observed_selcoefs_flat.reshape(-1, 1), predicted_selcoefs_flat)
    plt.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=plt.gca().transAxes, 
            verticalalignment='top')

    plt.tight_layout()
    plt.show()

# Permutation test - strategy 1
def permutation_test_type_1(
        A: np.ndarray,
        b: np.ndarray,
        obs_test_statistic: float, test_statistic: str = 'normalized_squared_error',
        num_permutations: int = 1000) -> Tuple[float, np.ndarray, float]:
    """
    Function to perform a permutation test on the fitness model.
    For the moment, only a normalized error squared is used as the metric.
    Strategy 1: take the vector of observations b min and max and sample uniformly from them.
    """

    # Get the min and max from the observation vector:
    b_min = np.min(b)
    b_max = np.max(b)

    # Get the length of the observation vector:
    b_len = len(b)

    # Initialize the vector of permutations
    permutations = np.zeros(num_permutations)

    for i in range(num_permutations):

        # Sample uniformly from min_b to max_b
        b_perm = np.random.uniform(b_min, b_max, b_len)

        # Solve the least squares problem
        solve_perm = solve_least_squares(A, b_perm)
        permutations[i] = solve_perm[test_statistic]

    # Calculate the p-value
    p_value = np.sum(permutations <= obs_test_statistic) / num_permutations

    # Count the number of permutation < 0:
    prob_neg = np.sum(permutations < 0) / num_permutations

    return p_value, permutations, prob_neg


def permutation_test_type_2(
        A: np.ndarray,
        m: np.ndarray,
        obs_test_statistic: float, test_statistic: str = 'normalized_squared_error',
        num_permutations: int = 1000) -> Tuple[float, np.ndarray, float]:
    """
    Function to perform a permutation test on the fitness model.
    For the moment, only a normalized error squared is used as the metric.
    Strategy 2: Take the matrix of observations m min and max and sample uniformly from them.
    """

    # Get the min and max from the matrix m:
    m_flat = m[m != 0]
    m_min = np.min(m_flat)
    m_max = np.max(m_flat) # or (-1) * m_min?

    # Get the number of b permuted elements:
    # This is number of elements in the matrix upper triangle.
    m_upper_len = int((m.shape[0] * (m.shape[0] - 1)) / 2)

    # Initialize the vector of permutations
    permutations = np.zeros(num_permutations)

    for i in range(num_permutations):

        # Sample uniformly from m_min to m_max
        b_perm = np.random.uniform(m_min, m_max, m_upper_len)

        # Solve the least squares problem
        solve_perm = solve_least_squares(A, b_perm)
        permutations[i] = solve_perm[test_statistic]

    # Calculate the p-value
    p_value = np.sum(permutations <= obs_test_statistic) / num_permutations

    # Count the number of permutation < 0:
    prob_neg = np.sum(permutations < 0) / num_permutations

    return p_value, permutations, prob_neg


def positive_permutation_wrapper(
    permutation_func: Callable[..., Tuple[float, np.ndarray, float]],
    n_positive: int = 1000,
    max_iterations: int = 10000,
    **kwargs: Any
) -> Tuple[List[float], int]:
    """
    Wrapper function to generate a specified number of positive permutation values.

    Args:
    permutation_func (Callable): The permutation test function to wrap.
    n_positive (int): Number of positive permutations to collect (default 1000).
    max_iterations (int): Maximum number of iterations to prevent infinite loops (default 10000).
    **kwargs: Additional keyword arguments to pass to the permutation function.

    Returns:
    Tuple[List[float], float, int, float]: A tuple containing the list of positive permutations,
                                           the calculated p-value, the total iterations performed,
                                           and the non-negative observed test statistic used.
    """
    obs_test_statistic = kwargs.pop('obs_test_statistic', None)
    if obs_test_statistic is None:
        raise ValueError("obs_test_statistic must be provided in kwargs")

    # Ensure non-negative observed test statistic
    obs_test_statistic = max(0, obs_test_statistic)
    
    positive_permutations = []
    iterations = 0

    while len(positive_permutations) < n_positive and iterations < max_iterations:
        _, permutation, _ = permutation_func(obs_test_statistic=obs_test_statistic, **kwargs)
        positive_values = permutation[permutation >= 0].tolist()
        positive_permutations.extend(positive_values)
        iterations += 1

        if iterations % 100 == 0:  # Print progress every 100 iterations
            print(f'Iteration {iterations}: {len(positive_permutations)} positive values collected')

    positive_permutations = positive_permutations[:n_positive]  # Trim to exact number if we've exceeded

    if iterations == max_iterations:
        print(f"Warning: Reached maximum iterations ({max_iterations}) before collecting {n_positive} positive values.")
    else:
        print(f"Collected {len(positive_permutations)} positive values in {iterations} iterations.")

    # Calculate p-value
    p_value = sum(perm <= obs_test_statistic for perm in positive_permutations) / len(positive_permutations)

    return positive_permutations, p_value, iterations, obs_test_statistic


def qq_distribution_plot(
        permuted_data: np.ndarray, 
        df: int,
        title: str = "Q-Q Distribution Plot") -> None:
    """
    Plot a Q-Q plot distribution plot comparing 
    the permuted data to a theoretical Chi-square distribution.

    Parameters:
    permuted_data (array-like): Array of permuted normalized squared errors.
    df (int): Degrees of freedom for the Chi-square distribution.
    title (str): Title for the combined plot.
    """

    # Sort the permuted data for the Q-Q plot
    sorted_permuted_data = np.sort(permuted_data)

    # Generate theoretical quantiles
    theoretical_quantiles = stats.chi2.ppf((np.arange(1, len(permuted_data) + 1) - 0.5) / len(permuted_data), df)
 
    # Create the scatter plot
    plt.figure(figsize=(10, 8))

    # Q-Q Plot
    plt.scatter(theoretical_quantiles, sorted_permuted_data, label='Permuted Data Quantiles', color='blue')
    plt.plot([theoretical_quantiles[0], theoretical_quantiles[-1]], [theoretical_quantiles[0], theoretical_quantiles[-1]], 'r--', label='Expected: Chi-square Quantiles')
    plt.xlabel('Theoretical Quantiles (Chi-square)', fontsize=14)
    plt.ylabel('Permuted Data Quantiles', fontsize=14)
    plt.legend()
    plt.title(title)
    plt.grid(False)
    plt.show()
    

def cdf_plot(
        permuted_data: np.ndarray, 
        df: int,
        title: str = "CDF Plot") -> None:
    """
    Plot a CDF distribution plot comparing 
    the permuted data to a theoretical Chi-square distribution.

    Parameters:
    permuted_data (array-like): Array of permuted normalized squared errors.
    df (int): Degrees of freedom for the Chi-square distribution.
    title (str): Title for the combined plot.
    """

    # Sort the permuted data for the CDF plot
    sorted_permuted_data = np.sort(permuted_data)
    cprob = np.array([i/len(sorted_permuted_data) for i in range(1,len(sorted_permuted_data)+1)])

    # Generate theoretical quantiles
    x = np.linspace(0, max(sorted_permuted_data), 200)

    # Create the scatter plot
    plt.figure(figsize=(10, 8))

    # cdf plot
    plt.plot(sorted_permuted_data, cprob, label=r'Observed: cumulative distribution of Normized Squared Error')
    plt.plot(x, stats.chi2.cdf(x, df=df), label=rf'Expected: $\chi^2 \ {{{df}}}\mathrm{{df}}$ ')
    plt.xlabel(rf'NSE, $\chi^2 \, {{{df}}}\mathrm{{df}}$', fontsize=14)
    plt.ylabel('Cumulative Probability', fontsize=14)
    plt.legend()
    plt.title(title)
    plt.grid(False)
    plt.show()


def histograma_plot(
        permuted_data: np.ndarray,
        obs_test_statistic: float,
        title: str = "Histograma Plot") -> None:
    """
    Plot a histogram distribution plot comparing 
    the permuted data to a theoretical Chi-square distribution.

    Parameters:
    permuted_data (array-like): Array of permuted normalized squared errors.
    obs_test_statistic (float): Observed test statistic.
    title (str): Title for the combined plot.
    """

    # Create the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(permuted_data, bins=20, density=True, alpha=0.6, color='b')
    plt.axvline(obs_test_statistic, color='r', linestyle='dashed', linewidth=2)

    # Add labels and title
    plt.xlabel('Normalized Squared Error', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.title(title)
    plt.grid(False)
    plt.show()


def fisher_method(
        p_values: List[float],
        df: int) -> Tuple[float, float]:
    """
    Perform Fisher's Methdo to calculate the chi-squared statistic and p-value.
    for a list of p-values.

    Parameters:
    p_values (list): List of p-values.

    Returns:
    Tuple of chi-squared statistic and p-value.
    """

    # calculate the chi-squared statistic
    chi_squared = -2 * np.sum(np.log(p_values))

    # calculate the p-value
    p_value = 1 - stats.chi2.cdf(chi_squared, df)

    return chi_squared, p_value


def main() -> None:
    """
    Main function to test the fitness_analyses function
    """

    # An example usage:

    # Input data
    codon_changes = [
        ('CCT>CCC', 0.129),
        ('CCA>CCT',-0.542),
        ('CCT>CCG',0.646),
        ('CCA>CCC',-0.155),
        ('CCG>CCC',0.087),
        ('CCA>CCG',-0.167)
    ]

    # Define the codons in order
    # Here you can define a different reference codon by setting it to index 0
    codons = ['CCA', 'CCC', 'CCG', 'CCT']

    # Create the m x m matrix of codon changes 2Ns'.
    m = create_codon_matrix(codon_changes, codons)

    # Build the design matrix and observation vector
    A = build_design_matrix(m.shape[0])

    # Build the observation vector
    b = build_observation_vector(m)

    # Solve the least squares problem
    solved = solve_least_squares(A, b)
    print(solved)

    # Create the matrix m x m of codon changes 2Ns'.
    # m_pred = construct_full_predicted_matrix(solved['predicted_selection_coefficients'])

    # # Plot selection coefficients observed vs predicted
    # scatterplot_selcoef_matrices(codons, m, m_pred)

    # Plot selection coefficients observed vs predicted
    scatterplot_selcoef_vectors(codons, b, solved['predicted_selection_coefficients'])
    
    # Permutation test type 1:
    # positive_values, p_value, total_iterations, used_obs_statistic = positive_permutation_wrapper(
    #     permutation_test_type_1,
    #     A=A,
    #     b=b,
    #     obs_test_statistic=solved['normalized_squared_error']
    # )

    # print(f"P-value: {p_value}, total iterations: {total_iterations}, used obs stat: {used_obs_statistic} for {len(positive_values)} permutations using type 1 test.")

    # Permutation test type 2:
    positive_values, p_value, total_iterations, used_obs_statistic = positive_permutation_wrapper(
        permutation_test_type_2,
        A=A,
        m=m,
        obs_test_statistic=solved['normalized_squared_error']
    )

    print(f"P-value: {p_value}, total iterations: {total_iterations}, used obs stat: {used_obs_statistic} for {len(positive_values)} permutations using type 2 test.")

    # CDF plot
    cdf_plot(positive_values, df=3, title="CDF Plot")

    # Histograma plot
    histograma_plot(positive_values, obs_test_statistic=solved['normalized_squared_error'], title="Histograma Plot")


if __name__ == '__main__':
    main()
