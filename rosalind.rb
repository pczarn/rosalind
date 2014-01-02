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

def iprb(dd, dr, rr)

end

def tree(n, edges)

end

def mprt(proteins)

end

def prob(dna, gc_contents)

end

def tran(dna1, dna2)

end

def protstruct(points)

end
