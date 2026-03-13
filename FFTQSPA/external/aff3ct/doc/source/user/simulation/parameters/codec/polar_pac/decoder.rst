.. _dec-polar_pac-decoder-parameters:

Polar |PAC| Decoder parameters
------------------------------

.. _dec-polar_pac-dec-type:

``--dec-type, -D``
""""""""""""""""""

   :Type: text
   :Allowed values: ``SCL``
   :Default: ``SCL``
   :Examples: ``--dec-type SCL``

|factory::Decoder::p+type,D|

Description of the allowed values:

+--------------+---------------------------------------------------------------+
| Value        | Description                                                   |
+==============+===============================================================+
| ``SCL``      | Select the |SCL| algorithm from :cite:`Yao2020`.              |
+--------------+---------------------------------------------------------------+

.. _dec-polar_pac-dec-implem:

``--dec-implem``
""""""""""""""""

   :Type: text
   :Allowed values: ``NAIVE``
   :Default: ``NAIVE``
   :Examples: ``--dec-implem NAIVE``

|factory::Decoder::p+implem|

Description of the allowed values:

+-----------+--------------------------+
| Value     | Description              |
+===========+==========================+
| ``NAIVE`` | |dec-implem_descr_naive| |
+-----------+--------------------------+

.. |dec-implem_descr_naive| replace:: Select the naive implementation which is
   typically slow.

.. _dec-polar_pac-dec-lists:

``--dec-lists, -L``
"""""""""""""""""""

   :Type: integer
   :Default: ``8``
   :Examples: ``--dec-lists 128``

|factory::Decoder_polar_PAC::p+lists,L|

References
""""""""""

.. bibliography:: references_dec.bib
   :labelprefix: Pacd-