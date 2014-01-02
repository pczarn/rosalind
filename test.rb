require "test/unit"

require_relative "rosalind"

class TestRosalind < Test::Unit::TestCase
   def test_dna_counting_nucleotides
      assert_equal([20, 12, 17, 21],
         dna("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))
   end

   def test_rna_transcribing_into_dna
      assert_equal("GAUGGAACUUGACUACGUAAAUU",
         rna("GATGGAACTTGACTACGTAAATT"))
   end

   def test_revc_completing_a_strand_of_dna
      assert_equal("ACCGGGTTTT",
         revc("AAAACCCGGT"))
   end

   def test_fib_rabbits_and_recurrence_relations
      assert_equal(19, fib(5, 3))
   end

   def test_gc_content_computing
      input = <<-INPUT
>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT
      INPUT
      id, gc_content = gc(input)
      assert_equal("Rosalind_0808", id)
      compare(60.919540, gc_content)
   end

   def test_hamm
      assert_equal(7, hamm("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"))
   end

   def test_iprb_mendels_first_law
      compare(0.78333, iprb(2, 2, 2))
   end

   def test_prot_translating_rna_into_protein
      assert_equal("MAMAPRTEINSTRING",
         prot("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
   end

   def test_subs_finding_a_motif_in_dna
      assert_equal([2, 4, 10],
         subs("GATATATGCATATACTT", "ATAT"))
   end

   def test_cons_consensus_and_profile
      input = <<-INPUT
>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT
      INPUT
      assert_equal(<<-OUT, cons(input))
ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6
      OUT
   end

   def test_fibd_mortal_fibonacci_rabbits
      assert_equal(4, fibd(6, 3))
   end

   def test_grph_overlap
      input = <<-INPUT
>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG
      INPUT
      assert_equal(<<-OUT, grph(input))
Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323
      OUT
   end

   def iev_calculating_expected_offspring
      compare(3.5, iev(1, 0, 0, 1, 0, 1))
   end

   def test_lcsm_finding_a_shared_motif
      assert_equal("AC", <<-IN)
>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA
      IN
   end

   def test_tree_completing
      input = [[1, 2], [2, 8], [4, 10], [5, 9], [6, 10], [7, 9]]
      assert_equal 3,
         tree(10, input)
   end

   def test_prob_gc_content_probability
      input = %w{0.129 0.287 0.423 0.476 0.641 0.742 0.783}.map(&:to_f)
      compare("-5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009",
         prob("ACGATACAA", input))
   end

   def test_tran_transitions_transversions
      input = %w(GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT
                 TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT)
      compare(1.21428571429,
         tran(*input))
   end

   # Suggested problems

   def test_distance_between_atoms
      input = %w{
         17.426 -32.764 65.278
         2.109 -38.295 41.517
         -7.758 -42.568 33.470
         -5.866 -3.238 13.986
         -4.720 -38.377 -4.862
      }.map(&:to_f).each_slice(3).to_a
      compare("28.806 41.738 63.601 73.767 13.430 45.283 46.879 43.932 38.680 39.891",
         protstruct(input))
   end

   def compare(expected, actual, precision=3)
      delta = 0.1**precision
      expected = expected.split.map(&:to_f) if expected.kind_of?(String)
      expected = Array(expected)
      actual = Array(actual)

      actual = actual.zip(expected).map do |act, exp|
         if act.kind_of?(Float) && exp.kind_of?(Float) && (exp - act).abs < delta
            exp
         else
            act
         end
      end

      assert_equal expected, actual
   end
end
