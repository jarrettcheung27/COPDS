.. _enc-polar_pac-encoder-parameters:

Polar |PAC| Encoder parameters
------------------------------

.. _enc-polar_pac-enc-type:

``--enc-type``
""""""""""""""

   :Type: text
   :Allowed values: ``POLAR_PAC`` ``AZCW`` ``COSET`` ``USER``
   :Default: ``POLAR_PAC``
   :Examples: ``--enc-type AZCW``

|factory::Encoder::p+type|

Description of the allowed values:

+---------------+------------------------+
| Value         | Description            |
+===============+========================+
| ``POLAR_PAC`` | |enc-type_descr_polar| |
+---------------+------------------------+
| ``AZCW``      | |enc-type_descr_azcw|  |
+---------------+------------------------+
| ``COSET``     | |enc-type_descr_coset| |
+---------------+------------------------+
| ``USER``      | |enc-type_descr_user|  |
+---------------+------------------------+

.. |enc-type_descr_polar| replace:: Select the standard Polar |PAC| encoder.
  The :ref:`enc-polar_pac-enc-conv` argument is required.
.. |enc-type_descr_azcw| replace:: See the common :ref:`enc-common-enc-type`
   parameter.
.. |enc-type_descr_coset| replace:: See the common :ref:`enc-common-enc-type`
   parameter.
.. |enc-type_descr_user| replace:: See the common :ref:`enc-common-enc-type`
   parameter.

.. note:: Only **non-systematic encoding** is supported as defined in
  :cite:`Arikan2019bis`.

.. _enc-polar_pac-enc-conv:

``--enc-conv``
""""""""""""""

   :Type: text
   :Default: ``NO``
   :Examples: ``--enc-conv 0o33``

|factory::Encoder_polar_PAC::p+conv|

.. warning:: It is required to specify this argument with the ``POLAR_PAC``
  encoder. 
.. note:: Currently, we only support the octal format for the impulse response of the convolutional transform. For instance, consider an impulse response represented as ``[1, 0, 1, 1, 1, 0, 1]``. When converting this to octal format, we will ignore the first ``1`` and convert the remaining vector, since ``AFF3CT`` will append the ``1`` at the start.

Thus, the vector ``[0, 1, 1, 1, 0, 1]`` will be converted to the octal value ``35``. Consequently, the input for ``AFF3CT`` will be ``0o35``.

.. _enc-polar_pac-enc-fb-gen-method:

``--enc-fb-gen-method``
"""""""""""""""""""""""

   :Type: text
   :Allowed values: ``RM`` ``FILE`` ``GA`` ``TV`` ``BEC`` ``5G``
   :Examples: ``--enc-fb-gen-method FILE``

|factory::Frozenbits_generator::p+gen-method|

Description of the allowed values:

+----------+-------------------------------------------------------------------+
| Value    | Description                                                       |
+==========+===================================================================+
| ``RM``   | This is the preferred method to generate short length |PAC| codes.|
+----------+-------------------------------------------------------------------+
| ``GA``   | Select the |GA| method from :cite:`Trifonov2012`.                 |
+----------+-------------------------------------------------------------------+
| ``TV``   | Select the |TV| method which is based on Density Evolution (|DE|) |
|          | approach from :cite:`Tal2013`, to use with the                    |
|          | :ref:`enc-polar_pac-enc-fb-awgn-path` parameter.                  |
+----------+-------------------------------------------------------------------+
| ``FILE`` | Read the best channels from an external file, to use with the     |
|          | :ref:`enc-polar_pac-enc-fb-awgn-path` parameter.                  |
+----------+-------------------------------------------------------------------+
| ``BEC``  | Generate frozen bits for the |BEC| channel from                   |
|          | :cite:`Arikan2009`.                                               |
+----------+-------------------------------------------------------------------+
| ``5G``   | Generate the frozen bits as described in the 5G standard          |
|          | :cite:`3GPP2017`.                                                 |
+----------+-------------------------------------------------------------------+

.. note:: The ``RM`` method is typically good for short block lengths and with
   PAC codes :cite:`Arikan2019bis`, :cite:`Yao2020`. This is the default
   generator on |PAC| codes.

.. note:: By default, when using the |GA| or the |TV| method, the frozen bits
   are optimized for each |SNR| point. To override this behavior you can use
   the :ref:`enc-polar_pac-enc-fb-noise` parameter.

.. note:: When using the ``FILE`` method, the frozen bits are always the same
   regardless of the |SNR| value.

.. note:: When using the ``BEC`` method, the frozen bits are optimized for each
   erasure probability.

.. note:: When using the ``5G`` method, the codeword size must be inferior
   to 1024.

.. _enc-polar_pac-enc-fb-awgn-path:

``--enc-fb-awgn-path``
""""""""""""""""""""""

   :Type: path
   :Rights: read only
   :Examples: ``--enc-fb-awgn-path example/path/to/the/right/place/``

|factory::Frozenbits_generator::p+awgn-path|

An |ASCII| file is expected, for instance, the following file describes the
most reliable channels optimized for a codeword of size :math:`N = 8` and for an
|AWGN| channel where the noise variance is :math:`\sigma = 0.435999`:

.. code-block:: bash

   8
   awgn
   0.435999
   7 6 5 3 4 2 1 0

Given the previous file, if we suppose a Polar code of size :math:`N = 8` with
:math:`K = 4` information bits, the frozen bits are at the ``0, 1, 2, 4``
positions in the codeword. The strategy is to freeze the less reliable channels.

.. warning:: The ``FILE`` frozen bits generator expects a file and not a
   directory.

.. warning:: The ``TV`` frozen bits generator expects a directory and not a
   file. |AFF3CT| comes with input configuration files, a part of those
   configuration files are a set of best channels pre-generated with the |TV|
   method (see ``conf/cde/awgn_polar_codes/TV/``).

.. _enc-polar_pac-enc-fb-dump-path:

``--enc-fb-dump-path``
""""""""""""""""""""""

   :Type: folder
   :Rights: write only
   :Examples: ``--enc-fb-dump-path example/path/to/the/right/place/``

|factory::Frozenbits_generator::p+dump-path|

.. note:: Works only for the ``GA`` and ``BEC`` frozen bits generation methods.

.. _enc-polar_pac-enc-fb-noise:

``--enc-fb-noise``
""""""""""""""""""

   :Type: real number
   :Examples: ``--enc-fb-noise 1.0``

|factory::Frozenbits_generator::p+noise|

Can be a gaussian noise variance :math:`\sigma` for |GA| and |TV| generation
methods, or an event probability for the |BEC| generation method. All the noise
points in the simulation will use the same frozen bits configuration.

References
""""""""""

.. bibliography:: references_enc.bib
   :labelprefix: Pace-
