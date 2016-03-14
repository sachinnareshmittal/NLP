#Assigment 2 Parag,Soham
########################################################################


###traing data
unigrams = {}
tagwords = {}

bigrams = {}
trigrams = {}


def bigram_count(sentence):
    for i in range(len(sentence)-1):
        u = sentence[i]
        v = sentence[i+1]
        if (u,v) in bigrams:
            bigrams[u,v] += 1
        else:
            bigrams[u,v] = 1

def trigram_count(sentence):
    for i in range(len(sentence)-2):
        u = sentence[i]
        v = sentence[i+1]
        w = sentence[i+2]
        if (u,v,w) in trigrams:
            trigrams[u,v,w] += 1
        else:
            trigrams[u,v,w] = 1


for line in open('Brown_tagged_train.txt','r').readlines():
    myline = line.split()
    sentence_tags = []
    for word in myline:
        tags = word.split('/')
        tag =  tags[-1]
        if tag in unigrams:
            unigrams[tag] += 1
        else:
            unigrams[tag] = 1
        word = ""
        for i in range(len(tags) - 1):
            word = word + tags[i] + '/'
        word = word[:-1]
        if (tag,word) in tagwords:
            tagwords[tag,word] += 1
        else:
            tagwords[tag,word] = 1
            
        sentence_tags.append(tag)
    sentence_tags = ['*'] + ['*'] + sentence_tags
    #print(sentence_tags)
    bigram_count(sentence_tags)
    trigram_count(sentence_tags)




def QML_s_given_u_v(s,u,v):
    if (u,v,s) in trigrams:
        num = trigrams[u,v,s]
        den = bigrams[u,v]
        f = (float)(num/den)
        return f
    else:
        return 0
    
def E_tag_to_word(tag,word):
    if (tag,word) in tagwords:
        num = tagwords[tag,word]
        den = unigrams[tag]
        f = (float)(num/den)
        return f
    else:
        return 0

def Q_s_given_u_v(s,u,v):
    return QML_s_given_u_v(s,u,v)
    


start_symbol = {}
end_symbol = {}
start_symbol['*'] = 'dummy'
end_symbol['.'] = 'dummy'


def veterbi(sentence):
    pi = {}
    pi[0,'*','*'] = 1
    bp = {}
    length = len(sentence)
    for k in range(1,length):
        possible_S_k_2 = start_symbol
        possible_S_k_1 = start_symbol
        possible_S_k = start_symbol
        if k == 1:
            possible_S_k_2 = start_symbol
            possible_S_k_1 = start_symbol
            possible_S_k = unigrams
        else:
            if k == 2:
                possible_S_k_2 = start_symbol
                possible_S_k_1 = unigrams
                possible_S_k = unigrams
            else:
                if k == length:
                    possible_S_k_2 = unigrams
                    possible_S_k_1 = unigrams
                    possible_S_k = end_symbol
                else:
                    possible_S_k_2 = unigrams
                    possible_S_k_1 = unigrams
                    possible_S_k = unigrams

        for v, value in possible_S_k.items():
            for u, value in possible_S_k_1.items():
                max_value = -1
                for w, value in possible_S_k_2.items():
                    now = pi[k-1,w,u] * Q_s_given_u_v(v,w,u) * E_tag_to_word(v,sentence[k-1])
                    if now > max_value:
                        max_value = now
                        bp[k,u,v] = w

                pi[k,u,v] = max_value
        #print(k)
        #print(length)

    tags = []
    maxi = -1
    ans = 1
    for i in range(0,length):
        tags.append('pp')
    #print(tags)
    for v, value in possible_S_k.items():
            for u, value in possible_S_k_1.items():
                now = (pi[length-1,u,v]* Q_s_given_u_v('.',u,v))
                if now > maxi:
                       ans = u,v
    tags[length - 1] = ans[0]
    tags[length - 2] = ans[1]
    kk = length -3
    while ( kk >= 0):
        #print(kk)
        one = tags[kk+1]
        two = tags[kk+2]
        tags[kk] = bp[kk+2,one ,two]
        kk = kk - 1
    return tags
        
            
    
                       
            

    



    




#testingaing data

for line in open('train.txt','r').readlines():
    myline = line.split()
    #I have to find a tag seq which maximises the probability
    tag_seq = veterbi(myline)
    print(tag_seq)
    
