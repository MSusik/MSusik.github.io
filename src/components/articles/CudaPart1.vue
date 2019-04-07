<template>
<div>
<h1>Implementing Smith-Waterman algorithm on CUDA</h1>
<h2>Part 1: C++ implementation</h2>
<p><em>07.04.2019</em></p>
<p>In this series of posts, I will show how to accelerate a heavy algorithm using a graphics processing unit
(GPU). The algorithm of my choice is Smith-Waterman. It suits me as its complexity is <img src="https://tex.s2cms.ru/svg/%5Cmathcal%7BO%7D(n%5E2)" alt="\mathcal{O}(n^2)" /> in both space and time. Firstly, I will present a C++ implementation. While it is possible to use CUDA from other programming languages (<a href="https://developer.nvidia.com/pycuda">PyCUDA</a>, <a href="http://www.jcuda.org">JCUDA</a>) I chose to use C++ for consistency with the GPU kernels. CUDA’s programming language is a superset of C, although there exists also a Fortran version. Note that you do not need to have an Nvidia graphics card in order to do computations on GPU - alternatively, <a href="https://en.wikipedia.org/wiki/OpenCL">OpenCL</a> can be successfully used on most of the platforms.</p>
<p>Having said that, let us quickly outline the algorithm.</p>
<p>Smith-Waterman algorithm aims to find the best local alignment of two protein or nucleic acid sequences. The user inputs the alignment reward <img src="https://tex.s2cms.ru/svg/s" alt="s" /> and gap penalty <img src="https://tex.s2cms.ru/svg/d" alt="d" />. Then, the algorithm computes local alignment scores defined as:</p>
<p align="center" style="text-align: center;"><img align="center" src="https://tex.s2cms.ru/svg/F(i%2Cj)%3D%5Cmax%5Cleft%5C%7B%0A%5Cbegin%7Barray%7D%7Blr%7D%0A0%5C%5C%0AF(i-1%2Cj)-d%2C%5C%5C%0AF(i%2C%20j-1)-d%2C%5C%5C%0AF(i-1%2Cj-1)%2B%20s(x_%7Bi%7D%2Cy_%7Bj%7D)%0A%5Cend%7Barray%7D%0A%5Cright%5C%7D" alt="F(i,j)=\max\left\{
\begin{array}{lr}
0\\
F(i-1,j)-d,\\
F(i, j-1)-d,\\
F(i-1,j-1)+ s(x_{i},y_{j})
\end{array}
\right\}" /></p>
<p>where <img src="https://tex.s2cms.ru/svg/i" alt="i" /> and <img src="https://tex.s2cms.ru/svg/j" alt="j" /> are the indices of sequences. Once the local alignment scores are computed, in order to find the best local alignment, the algorithm starts from the maximal element and &quot;backtracks&quot; the path by traversing recursively to the source of the score. A more detailed explanation can be found on <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Wikipedia</a>.</p>
<h3>Class definition</h3>
<p>In order to organize my code well, I will define a <code>SmithWatermanAlgorithm</code> class that will be responsible for running the algorithm:</p>
<pre><code language="cpp">class SmithWatermanAlgorithm {
    private:
        std::vector&lt;int&gt; scoring_matrix;
        std::vector&lt;unsigned char&gt; backtracking_matrix;
        const std::string&amp; first, other;
        int scoring_size(const std::string&amp;, const std::string&amp;);
        int backtracking_size(const std::string&amp;, const std::string&amp;);
        std::tuple&lt;int, int&gt; max_indices();
        int row_size;
        int gap_penalty;
    public:
        SmithWatermanAlgorithm (const std::string&amp;, const std::string&amp;, int);
        void fill_matrices(const substitution_matrix_t&amp;);
        void print_scoring_matrix();
        void print_backtracking_matrix();
        std::tuple&lt;std::string, std::string&gt; backtrack();
};
</code></pre>
<p>The constructor takes two sequences and the gap penalty as arguments. I went for a 2D array represented as a 1D array due to simplicity in case of the following GPU version. If you plan to develop only for CPU, it is easier to read the code that works on proper 2D arrays and follows modern C++, eg. the <a href="https://en.cppreference.com/w/cpp/language/range-for">range-based for</a> construct.</p>
<pre><code language="cpp">
SmithWatermanAlgorithm::SmithWatermanAlgorithm (
        const std::string&amp; first,
        const std::string&amp; other,
        int gap_penalty
):
    scoring_matrix(this-&gt;scoring_size(first, other), 0),
    backtracking_matrix(this-&gt;backtracking_size(first, other), 0),
    first(first),
    other(other),
    row_size(first.size() + 1),
    gap_penalty(gap_penalty)
{}

int SmithWatermanAlgorithm::scoring_size(
        const std::string&amp; first,
        const std::string&amp; other
) {
    # One cell padding from left and top
    return (first.size() + 1) * (other.size() + 1);
}

int SmithWatermanAlgorithm::backtracking_size(
        const std::string&amp; first,
        const std::string&amp; other
) {
    return first.size() * other.size();
}
</code></pre>
<h2>Scoring matrix</h2>
<p>The algorithm starts with the creation of a scoring matrix. I will also create a separate backtracking matrix that will store the information about the gaps used. It will be useful in the CUDA scenario. The backtracking algorithm will be run on CPU (as it is linear) and only the backtracking matrix will be copied back from the GPU in order to limit the copy size. Often, in the case of GPU applications, the cost associated with the copying procedures is dominating the computation cost.</p>
<pre><code language="cpp">void SmithWatermanAlgorithm::fill_matrices(const substitution_matrix_t&amp; substitution_matrix) {
    for ( int i = 0; i &lt; this-&gt;row_size - 1; i ++ )
    {
        for ( int j = 0; j &lt; this-&gt;other.size(); j ++ )
        {
            int diagonal_index = j * this-&gt;row_size + i;
            int top_index = diagonal_index + this-&gt;row_size;
            int left_index = diagonal_index + 1;
            int current = top_index + 1;
            int backtracking_index = diagonal_index - j;

            int diagonal_value = maximum(
                scoring_matrix[diagonal_index] + substitution_matrix.at({
                    this-&gt;first[i],
                    this-&gt;other[j]
                }),
                0
            );
            int left_value = maximum(scoring_matrix[left_index] - gap_penalty, 0);
            int top_value = maximum(scoring_matrix[top_index] - gap_penalty, 0);

            int max = maximum(diagonal_value, left_value, top_value);

            this-&gt;scoring_matrix[current] = max;
            if(max &gt; 0)
            {
                if (diagonal_value == max) {
                    this-&gt;backtracking_matrix[backtracking_index] = 'D';
                }
                else if (top_value == max) {
                    this-&gt;backtracking_matrix[backtracking_index] = 'T';
                }
                else if (left_value == max) {
                    this-&gt;backtracking_matrix[backtracking_index] = 'L';
                }
            }
        }
    }
}
</code></pre>
<p>In the backtracking matrix, the only information stored is the direction. It is enough for us to reconstruct the alignment later.</p>
<h3>Backtracking</h3>
<p>Once the backtracking matrix has been computed, we can run the procedure that will find the best local alignment. To do that, we find the maximal element:</p>
<pre><code language="cpp">std::tuple&lt;int, int&gt; SmithWatermanAlgorithm::max_indices() {
    int max_i = 0;
    int max_j = 0;
    int max_value = 0;
    for ( int i = 1; i &lt;= this-&gt;row_size - 1; i ++ )
    {
        for ( int j = 1; j &lt;= this-&gt;other.size(); j ++ )
        {
            int current_value = this-&gt;scoring_matrix[i + j * this-&gt;row_size];
            if (current_value &gt; max_value) {
                max_value = current_value;
                max_i = i - 1;
                max_j = j - 1;
            }
        }
    }

    return std::make_tuple(max_i, max_j);
}
</code></pre>
<p>and starting from it we find the path of the best local alignment.</p>
<pre><code language="cpp">std::tuple&lt;std::string, std::string&gt; SmithWatermanAlgorithm::backtrack() {
    auto indices = this-&gt;max_indices();

    int current_i = std::get&lt;0&gt;(indices);
    int current_j = std::get&lt;1&gt;(indices);

    std::string first_string, second_string;

    char value;
    do {
        value = this-&gt;backtracking_matrix[current_i + current_j * (this-&gt;row_size - 1)];
        switch(value) {
            case 'D':
                first_string += this-&gt;first[current_i];
                second_string += this-&gt;other[current_j];
                current_i--;
                current_j--;
                break;
            case 'L':
                first_string += '-';
                second_string += this-&gt;other[current_j];
                current_j--;
                break;
            case 'T':
                first_string += this-&gt;first[current_i];
                second_string += '-';
                current_i--;
                break;
        }
    }
    while (value != 0);

    std::reverse(first_string.begin(), first_string.end());
    std::reverse(second_string.begin(), second_string.end());

    return std::make_tuple(first_string, second_string);
}
</code></pre>
<p>While traversing the matrix we build the pair of strings that represents the alignment. Obviously, after the algorithm, we have to reverse them.</p>
<h3>Example</h3>
<p>I will start by showing the algorithm on a small example. The alignment reward was 3 for all matching pairs, the gap penalty was 2.</p>
<pre><code language="cpp">&gt;Test1
CGGGTATCCAA
&gt;Test2
CCCTAGGTCCCA
</code></pre>
<p>I added some debugging infrastructure in order to check what the algorithm does. For the convenience of the readers, I present the results in the form of clear tables. The output of the algorithm is:</p>
<pre><code language="cpp">SCORING MATRIX
</code></pre>
<table>
<thead>
<tr>
<th></th>
<th>|</th>
<th style="text-align:center">C</th>
<th style="text-align:center">C</th>
<th style="text-align:center">C</th>
<th style="text-align:center">T</th>
<th style="text-align:center">A</th>
<th style="text-align:center">G</th>
<th style="text-align:center">G</th>
<th style="text-align:center">T</th>
<th style="text-align:center">C</th>
<th style="text-align:center">C</th>
<th style="text-align:center">C</th>
<th style="text-align:center">A</th>
</tr>
</thead>
<tbody>
<tr>
<td>—</td>
<td>+</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
<td style="text-align:center">—</td>
</tr>
<tr>
<td><strong>C</strong></td>
<td>|</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
</tr>
<tr>
<td><strong>G</strong></td>
<td>|</td>
<td style="text-align:center">1</td>
<td style="text-align:center">1</td>
<td style="text-align:center">1</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">1</td>
<td style="text-align:center">1</td>
<td style="text-align:center">1</td>
<td style="text-align:center">0</td>
</tr>
<tr>
<td><strong>G</strong></td>
<td>|</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
<td style="text-align:center">2</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
</tr>
<tr>
<td><strong>G</strong></td>
<td>|</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
<td style="text-align:center">2</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
</tr>
<tr>
<td><strong>T</strong></td>
<td>|</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">1</td>
<td style="text-align:center">4</td>
<td style="text-align:center">9</td>
<td style="text-align:center">7</td>
<td style="text-align:center">5</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
</tr>
<tr>
<td><strong>A</strong></td>
<td>|</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">1</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
<td style="text-align:center">2</td>
<td style="text-align:center">7</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
<td style="text-align:center">2</td>
<td style="text-align:center">6</td>
</tr>
<tr>
<td><strong>T</strong></td>
<td>|</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">4</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">5</td>
<td style="text-align:center">4</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">4</td>
</tr>
<tr>
<td><strong>C</strong></td>
<td>|</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">2</td>
<td style="text-align:center">1</td>
<td style="text-align:center">0</td>
<td style="text-align:center">3</td>
<td style="text-align:center">8</td>
<td style="text-align:center">7</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
</tr>
<tr>
<td><strong>C</strong></td>
<td>|</td>
<td style="text-align:center">3</td>
<td style="text-align:center">6</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
<td style="text-align:center">2</td>
<td style="text-align:center">0</td>
<td style="text-align:center">0</td>
<td style="text-align:center">1</td>
<td style="text-align:center">6</td>
<td style="text-align:center">11</td>
<td style="text-align:center">10</td>
<td style="text-align:center">8</td>
</tr>
<tr>
<td><strong>A</strong></td>
<td>|</td>
<td style="text-align:center">1</td>
<td style="text-align:center">4</td>
<td style="text-align:center">4</td>
<td style="text-align:center">3</td>
<td style="text-align:center">7</td>
<td style="text-align:center">5</td>
<td style="text-align:center">3</td>
<td style="text-align:center">1</td>
<td style="text-align:center">4</td>
<td style="text-align:center">9</td>
<td style="text-align:center">8</td>
<td style="text-align:center"><strong><strong><em>13</em></strong></strong></td>
</tr>
<tr>
<td><strong>A</strong></td>
<td>|</td>
<td style="text-align:center">0</td>
<td style="text-align:center">2</td>
<td style="text-align:center">2</td>
<td style="text-align:center">1</td>
<td style="text-align:center">6</td>
<td style="text-align:center">4</td>
<td style="text-align:center">2</td>
<td style="text-align:center">0</td>
<td style="text-align:center">2</td>
<td style="text-align:center">7</td>
<td style="text-align:center">6</td>
<td style="text-align:center">11</td>
</tr>
</tbody>
</table>
<pre><code language="cpp">BACKTRACKING MATRIX
</code></pre>
<table>
<thead>
<tr>
<th></th>
<th>|</th>
<th>C</th>
<th>C</th>
<th>C</th>
<th>T</th>
<th>A</th>
<th>G</th>
<th>G</th>
<th>T</th>
<th>C</th>
<th>C</th>
<th>C</th>
<th>A</th>
</tr>
</thead>
<tbody>
<tr>
<td>—</td>
<td>+</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
<td>—</td>
</tr>
<tr>
<td><strong>C</strong></td>
<td>|</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
</tr>
<tr>
<td><strong>G</strong></td>
<td>|</td>
<td>↓</td>
<td>↓</td>
<td>↓</td>
<td>⠀</td>
<td>⠀</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
<td>↓</td>
<td>↓</td>
<td>↓</td>
<td>⠀</td>
</tr>
<tr>
<td><strong>G</strong></td>
<td>|</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
<td>→</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
</tr>
<tr>
<td><strong>G</strong></td>
<td>|</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
<td>→</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
</tr>
<tr>
<td><strong>T</strong></td>
<td>|</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>↘</td>
<td>→</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>→</td>
<td>→</td>
<td>→</td>
<td>→</td>
</tr>
<tr>
<td><strong>A</strong></td>
<td>|</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>↓</td>
<td>↘</td>
<td>→</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
</tr>
<tr>
<td><strong>T</strong></td>
<td>|</td>
<td>⠀</td>
<td>⠀</td>
<td>⠀</td>
<td>↘</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↓</td>
</tr>
<tr>
<td><strong>C</strong></td>
<td>|</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>⠀</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
</tr>
<tr>
<td><strong>C</strong></td>
<td>|</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
<td>→</td>
<td>⠀</td>
<td>⠀</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
</tr>
<tr>
<td><strong>A</strong></td>
<td>|</td>
<td>↓</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
<td>→</td>
<td>→</td>
<td>→</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
</tr>
<tr>
<td><strong>A</strong></td>
<td>|</td>
<td>⠀</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>↘</td>
<td>⠀</td>
<td>↓</td>
<td>↓</td>
<td>↘</td>
<td>↘</td>
</tr>
</tbody>
</table>
<pre><code language="cpp">BEST ALIGNMENT
GGTATCCA
GGT-CCCA
</code></pre>
<p>The scoring matrix shows the absolute value of the alignment score. The backtracking matrix shows how these scores were achieved. Finally, the best local alignment is printed.</p>
<h3>Benchmarking</h3>
<p>In this series of posts, I will be using my personal notebook. The related hardware is:</p>
<ul>
<li>Intel® Core™ i7-4720HQ CPU @ 2.60GHz</li>
<li>GeForce GTX965M</li>
</ul>
<p>I used <code>-O3</code> and <code>-ffast-math</code> flags during compilation.</p>
<p>For the purposes of the benchmarking, I will use example DNA sequences from <a href="http://genome.crg.es/datasets/genomics96/#SEQS">http://genome.crg.es/datasets/genomics96/#SEQS</a>. The alignment of two longer sequences (ACU08131 and AGGGLINE) took averagely 1,493 seconds.</p>
<p>We will see if it is possible to beat that time using a graphics card.</p>
<h3>Remarks</h3>
<p>In the world of HPC (high-performance computing) there exist many different ways of code parallelization. Some of them include:</p>
<ul>
<li>Multithreading</li>
<li>Multiprocessing</li>
<li>Distributed computing</li>
<li>Asynchronosity (though this one helps mostly with I/O active waiting and therefore is more suited to applications like web-development).</li>
</ul>
<p>In order to find out if it is worth to implement the algorithm on a CUDA device, one would have to create an optimized CPU version firstly. And only if the results are unsatisfactory, and proper profiling has been performed, I would suggest exploring other options, such as device-dependant CUDA.</p>
<p>In the case of the Smith-Waterman algorithm, if one plans to run many alignments, the simplest form of parallelization is to run each pair alignment in a separate process. This does not require any code change. Note that there are more efficient alignment algorithms like <a href="http://globin.cse.psu.edu/courses/spring2000/DP.pdf">Myers-Miller’s</a> and you might prefer to use them instead.</p>
</div>
</template>

<script>

export default {
  mounted () {
    document.querySelectorAll('pre code').forEach((block) => {
      this.$hljs.highlightBlock(block)
    })
  },
  name: 'CudaPart1',
  data () {
    return {
      name: this.$route.params.id
    }
  }
}
</script>
