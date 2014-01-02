require "test/unit"

require_relative "rosalind"

class TestRosalind < Test::Unit::TestCase
   def test_completing_a_tree
      input = [[1, 2], [2, 8], [4, 10], [5, 9], [6, 10], [7, 9]]
      assert_equal 3,
         tree(10, input)
   end

   def test_gc_content_probability
      input = %w{0.129 0.287 0.423 0.476 0.641 0.742 0.783}.map(&:to_f)
      compare("-5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009",
         prob("ACGATACAA", input))
   end

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

   def test_transitions_transversions
      input = %w(GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT
                 TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT)
      compare(1.21428571429,
         tran(*input))
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
