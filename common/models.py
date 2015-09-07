from django.db import models
from Bio import Entrez, Medline
from string import Template

import urllib.request,json

import logging


class WebResource(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=200, default='')
    url = models.TextField()
    #url should be a string template, so it can be automaticaly filled with index in proper place
    #https://docs.python.org/3.4/library/string.html#string.Template
    #Example: 'http://www.ncbi.nlm.nih.gov/pubmed/$index'

    def __str__(self):
        return self.url
    
    class Meta():
        db_table = 'web_resource'

class WebLink(models.Model):
    web_resource = models.ForeignKey('WebResource')
    index = models.TextField()

    #And now generating the working url is just a piece of cake
    def __str__(self):
        return Template(str(self.web_resource)).substitute(index=self.index)
    
    class Meta():
        db_table = 'web_link'


class Publication(models.Model):
    web_link = models.ForeignKey('common.WebLink', null=True)
    journal = models.ForeignKey('PublicationJournal')
    title = models.TextField()
    authors = models.TextField()
    year = models.IntegerField()
    reference = models.TextField()

    def __str__(self):
        return "{!s} ({!s})".format(self.journal, self.year)

    class Meta():
        db_table = 'publication'

    def update_from_pubmed_data(self, index=None):
        logger = logging.getLogger('build')
        try:
            Entrez.email = 'info@gpcrdb.org'
            if index:
                handle = Entrez.efetch(
                    db="pubmed", 
                    id=str(index), 
                    rettype="medline", 
                    retmode="text"
                    )
            else:
                handle = Entrez.efetch(
                    db="pubmed", 
                    id=str(self.web_link.index),
                    rettype="medline",
                    retmode="text"
                    )
        except Exception as e:
            logger.error("Failed 1x to retrieve data for pubmed id: {} - trying again".format(index))
            try:
                Entrez.email = 'info@gpcrdb.org'
                if index:
                    handle = Entrez.efetch(
                        db="pubmed", 
                        id=str(index), 
                        rettype="medline", 
                        retmode="text"
                        )
                else:
                    handle = Entrez.efetch(
                        db="pubmed", 
                        id=str(self.web_link.index),
                        rettype="medline",
                        retmode="text"
                        )
            except Exception as e:
                logger.error("Failed 2x to retrieve data for pubmed id: {}".format(index))
                return
        try:
            record = Medline.read(handle)
            self.title = record['TI']
            self.authors = record['AU']
            self.year = record['DA'][:4]
            record['JT'] = record['JT'][:29] #not allowed longer by model. FIXME
            record['TA'] = record['TA'][:29] #not allowed longer by model. FIXME
            try:
                self.journal = PublicationJournal.objects.get(slug=record['TA'], name=record['JT'])
            except PublicationJournal.DoesNotExist:
                j = PublicationJournal(slug=record['TA'], name=record['JT'])
                j.save()
                self.journal = j
            
            self.reference = ""
            if 'VI' in record:
                self.reference += record['VI']
            if 'IP' in record:
                self.reference += "(" + record['IP'] + ")"
            if 'PG' in record:
                self.reference += ":" + record['PG']
        except Exception as msg:
            logger.error("Publication update on pubmed error! Pubmed: "+index+" error:" + str(msg) )

    def update_from_doi(self, doi):
        logger = logging.getLogger('build')
        #Pubmed by doi
        try:
            url = "http://search.crossref.org/dois?q={}".format(doi)
            pub = json.loads(urllib.request.urlopen(url).read().decode('ascii', "ignore"))
            self.title = pub[0]['title'].strip()
            self.year = pub[0]['year'].strip()
            full_cit = [x.strip() for x in reversed(pub[0]['fullCitation'].split(','))]
            tmp_authors = []   
            p = re.compile(r'\d\d\d\d') #FIXME SOMETHING ERRORS HERE
            item = full_cit.pop()  
            while item.strip()[0] != "'" and p.match(item.strip()) is None:
                tmp_authors.append(item.strip())
                item = full_cit.pop()
            self.authors = ','.join(tmp_authors)      
            pages = ''
            vol = ''
            issue = ''    
            for i in full_cit:
                if '<i>' in i:
                    try:
                        self.journal = PublicationJournal.objects.get(name=i.replace('<i>','').replace('</i>',''))
                    except PublicationJournal.DoesNotExist:
                        j = PublicationJournal(name=i.replace('<i>','').replace('</i>',''))
                        j.save()
                        self.journal = j
                if i.replace("'",'').strip() == title:
                    continue
                if i.startswith('pp.'):
                    pages = i.replace('pp. ','')
                if i.startswith('no.'):
                    issue = '({})'.format(i.replace('no. ',''))
                if i.startswith('vol.'):
                    vol = i.replace('vol. ','')
            self.reference = "{}{}:{}".format(vol,issue,pages)    
        except Exception as msg:
            #Crossref doi search
            logger.error("Publication update on doi error - trying pubmed! DOI: "+index+" error:" + str(msg) )
            try:
                Entrez.email = 'info@gpcrdb.org'
                record = Entrez.read(Entrez.esearch(
                    db='pubmed',
                    retmax=1,
                    term=doi
                    ))
                self.update_from_pubmed_data(record['IdList'][0])

            except Exception as e:
                logger.error("Publication update on pubmed error! DOI: "+index+" error:" + str(e) )
                pass


class PublicationJournal(models.Model):
    slug = models.CharField(max_length=30, null=True)
    name = models.TextField()

    def __str__(self):
        return "{!s} ({!s})".format(self.name, self.slug)

    class Meta():
        db_table = 'publication_journal'