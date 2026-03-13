.. _dec-rep-decoder-parameters:

Repetition Decoder parameters
-----------------------------

.. _dec-rep-dec-type:

``--dec-type, -D``
""""""""""""""""""

   :Type: text
   :Allowed values: ``REPETITION`` ``REPETITION_VOTE`` ``CHASE`` ``ML``
   :Default: ``REPETITION``
   :Examples: ``--dec-type CHASE``

|factory::Decoder::p+type,D|

Description of the allowed values:

+---------------------+-----------------------------+
| Value               | Description                 |
+=====================+=============================+
| ``REPETITION``      | |dec-type_descr_repetition| |
+---------------------+-----------------------------+
| ``REPETITION_VOTE`` | |dec-type_descr_rep_vote|   |
+---------------------+-----------------------------+
| ``CHASE``           | |dec-type_descr_chase|      |
+---------------------+-----------------------------+
| ``ML``              | |dec-type_descr_ml|         |
+---------------------+-----------------------------+

.. |dec-type_descr_repetition| replace:: Select the repetition decoder (computes
   the mean probability per information bit)
.. |dec-type_descr_rep_vote| replace:: Select the repetition decoder (performs a
   majority vote on hard decided inputs).
.. |dec-type_descr_chase| replace:: See the common :ref:`dec-common-dec-type`
   parameter.
.. |dec-type_descr_ml| replace:: See the common :ref:`dec-common-dec-type`
   parameter.

.. _dec-rep-dec-implem:

``--dec-implem``
""""""""""""""""

   :Type: text
   :Allowed values: ``STD`` ``FAST``
   :Default: ``STD``
   :Examples: ``--dec-implem FAST``

|factory::Decoder::p+implem|

Description of the allowed values:

+-----------+--------------------------+
| Value     | Description              |
+===========+==========================+
| ``STD``   | |dec-implem_descr_std|   |
+-----------+--------------------------+
| ``FAST``  | |dec-implem_descr_fast|  |
+-----------+--------------------------+

.. |dec-implem_descr_std| replace:: Select the |STD| implementation.
.. |dec-implem_descr_fast| replace:: Select the fast implementation, much more
   faster without the :ref:`enc-rep-enc-no-buff` parameter.