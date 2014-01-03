require "httparty"

class FASTA < Hash
   def initialize(input)
      input.each_line.chunk do |line|
         line.chomp!
         line.sub!(/\A>\s*/, '') ? true : false
      end.each_slice(2) {|((_, (id)), (_, seq))|
         self[id] = seq.join('')
      }
   end

   def output

   end
end

class UniProt
   include HTTParty
   base_uri 'uniprot.org'

   def protein(id)
     self.class.get("/uniprot/#{ id }.txt")
   end

   def protein_fasta(id)
     self.class.get("/uniprot/#{ id }.fasta")
   end

   def post(text)
     options = { :body => {:status => text}, :basic_auth => @auth }
     self.class.post('/statuses/update.json', options)
   end
end
