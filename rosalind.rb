class Integer
  # binomial coefficient: n C k
  def choose(k)
    # n!/(n-k)!
    pTop = (self-k+1 .. self).inject(1, &:*)
    # k!
    pBottom = (2 .. k).inject(1, &:*)
    pTop / pBottom
  end
end

def dna(s)

end

def rna(s)

end

def revc(s)

end

def fib(n, k)

end

def gc(fasta)

end

def hamm(s, t)

end

def iprb(k, m, n)

end

def prot(s)

end

def subs(s, t)

end

def cons(fasta)

end

def fibd(n, m)

end

def grph(fasta)

end

def iev(*population)

end

def lcsm(fasta)

end

def tree(n, edges)

end

def mprt(proteins)

end

def prob(string, gc_contents)

end

def tran(str1, str2)

end

def protstruct(points)

end
