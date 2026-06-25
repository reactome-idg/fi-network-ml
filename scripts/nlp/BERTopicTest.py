from bertopic import BERTopic
from sklearn.datasets import fetch_20newsgroups
from random import  sample

docs = fetch_20newsgroups(subset='test', remove=('headers', 'footers', 'quotes'))['data']

topic_model = BERTopic()
docs = sample(docs, 100)

for i in range(len(docs)):
    print(i)
    print(docs[i])

#
# topics, probs = topic_model.fit_transform(docs)
#
# print(topic_model.get_topic_info())
#
# print(topic_model.get_topic(0))
